import pandas as pd
import numpy as np
from urllib.request import Request, urlopen
from bs4 import BeautifulSoup
import os
import wget
import subprocess
import timeit
import time

'''
traversal function:
    iterate through Refseq online driectories, grabbing the correct links and current pages being displayed
'''
def traversal(taxa, spcies, ref_or_gen):
    # Get all links (which are the taxonomic items)
    input_search = spcies.replace(' ', '_')
    
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/{ref_or_gen}/" + taxa + '/' + input_search + '/latest_assembly_versions'

    req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    html = urlopen(req).read()
    soup = BeautifulSoup(html, 'lxml')

    all_links = soup.find_all("a")
    
    dir_of_interest = ''
    for i in str(all_links[1])[9:]:
        if i == '/':
            break
        dir_of_interest = dir_of_interest + i

    dir_of_interest = dir_of_interest.replace(' ', '_')
    # Get all links (which are the taxonomic items)
    new_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/{ref_or_gen}/" + taxa + '/' + input_search + '/latest_assembly_versions/' + dir_of_interest
    req = Request(new_url, headers={'User-Agent': 'Mozilla/5.0'})
    html = urlopen(req).read()
    soup = BeautifulSoup(html, 'lxml')

    further_links = soup.find_all("a")
    
    url = new_url + '/' + dir_of_interest + '_protein.faa.gz'

    file_n = dir_of_interest
    file_faa = dir_of_interest + '_protein.faa'
    # Check if file exists, and if not, download and/or unzip it using gzip
    if os.path.exists('./absense_fasta_downloads/' + file_faa):
        print(file_n + ' is already unzipped in your directory!')
        return spcies, file_n
    elif os.path.exists('./absense_fasta_downloads/' + dir_of_interest + '_protein.faa.gz'):
        print('Already downloaded ' + dir_of_interest + '_protein.faa.gz. Unzipping now.')
        zip_dir = dir_of_interest + '_protein.faa.gz'
        !gunzip ./absense_fasta_downloads/$zip_dir
        return spcies, file_n
    else:
        try:
            wget.download(url, './absense_fasta_downloads')
            print('Successfully downloaded ' + dir_of_interest + '_protein.faa.gz. Unzipping now.')
            zip_dir = dir_of_interest + '_protein.faa.gz'
            !gunzip ./absense_fasta_downloads/$zip_dir
            print('Successfully grabbed annotated FASTA file called ' + zip_dir + ' for species ' + input_search)
            return spcies, file_n
        except:
            return 0, 0
       
'''
create_hash function:
    loop through the online taxonomic groups and create a hash table for faster accession in future iterations
    allows the program to correctly traverse the online directories based on only the inputted species name
'''
def create_hash(ref_or_gen):
    # Iterate through online directories and construct a hash table (dictionary) from species to groupings
    # We can ignore archaea, bacteria, and viruses because abSENSE does not work correctly with their genomes
    taxa_groups = ['fungi', 'invertebrate', 'plant', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other']

    species_dict = {}
    
    if ref_or_gen == 'REF':
        ins = 'Refseq'
    else:
        ins = 'Genbank'
    
    for k in taxa_groups:
        # Get all links (which are the taxonomic items within each subgroup)
        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/{ins.lower()}/{k}/"
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        html = urlopen(req).read()
        soup = BeautifulSoup(html, 'lxml')
        
        all_links = soup.find_all("a")
        
        # Loop through and create a mapping from each species to its group
        for i in all_links:
            val = str(i).strip('<a href="')[0:-5].split('/">').pop().replace('_', ' ')
            if 'annotation' not in val and 'Parent' not in val and 'summary' not in val:
                species_dict[val] = k
                
    with open(f"./absense_fasta_downloads/SPECIES_TAXA_HASH_{ref_or_gen}.txt","w") as f:
        for i, j in species_dict.items():
            f.write(i + ', ' + j + '\n')
        f.close()
   
'''
open_hash function:
    loop through the locally downloaded HASH text file and transform it to a dictionary
'''
def open_hash(inp_file):
    # Iterate through SPECIES_TAXA_HASH.txt and create Python dictionary object
    d = {}
    with open(inp_file) as f:
        for line in f:
            key, value = line.strip().split(",")
            d[key] = value.strip(' ')
    return d

'''
get_link function:
    from the large assembly dataframe, grab the relevant accession URLs to grab the protein FASTA files
'''
def get_link(data_assembly, input_search):

    # Assure that the species with full genome representations are chosen, and get proper link
    spec_data = data_assembly[data_assembly['organism_name'].fillna('na').str.contains(input_search)]
    
    # Add condition for if species is not in database
    if len(spec_data) == 0:
        return 0, 0
    fin_data = spec_data[spec_data['genome_rep'] == 'Full']
    
    # Add condition for if genome is not complete
    if len(fin_data) == 0:
        return 0, 0

    if sum(fin_data['refseq_category'] == 'reference genome') > 0:
        dir_data = fin_data[fin_data['refseq_category'] == 'reference genome']
    elif sum(fin_data['refseq_category'] == 'representative genome') > 0:
        dir_data = fin_data[fin_data['refseq_category'] == 'representative genome']
    elif sum(fin_data['relation_to_type_material'].isna()) != len(list(fin_data['relation_to_type_material'])):
        dir_data = fin_data[fin_data['relation_to_type_material'].notna()]
    else:
        dtes = list(pd.to_datetime(list(fin_data['seq_rel_date'])))
        dir_data = fin_data[pd.to_datetime(fin_data['seq_rel_date']) == max(dtes)]
    
    ftp = dir_data['ftp_path'].iloc[0]
    
    if 'vir' in dir_data['organism_name'].iloc[0]:
        return 0, 0
    
    # Construct link for protein annotation
    acc_name = dir_data['# assembly_accession'].iloc[0] + '_' + dir_data['asm_name'].iloc[0]
    lnk = 'https' + ftp[3:] + '/' + acc_name.replace(' ', '_') + '_protein.faa.gz'
    return lnk, acc_name

'''
alt_search function:
    if species protein annotation was not found in the assembly dataframe, manually traverse the online directories
    in search of the recent assembly download
'''
def alt_search(input_search, ref_or_gen):
    # Check if Refseq or Genbank option
    if ref_or_gen == 'REF':
        ins = 'Refseq'
    else:
        ins = 'Genbank'
    
    if os.path.exists(f'./absense_fasta_downloads/SPECIES_TAXA_HASH_{ref_or_gen}.txt'):
        print(f'{ins} species hash mappings already created, traversing online directories now...')
    else:
        print(f'Creating species taxa hash mappings for {ins}...')
        create_hash(ref_or_gen)
      
    # Iterate through hash table and use it for a file traversal to check if annotation exists
    hash_tbl = open_hash(f"./absense_fasta_downloads/SPECIES_TAXA_HASH_{ref_or_gen}.txt")
    try:
        taxa_cat = hash_tbl[input_search]
        res_spc, res_acc = traversal(taxa_cat, input_search, ins.lower())
        if res_acc != 0:
            return res_spc, res_acc
        else:
            print(f'Alternative {ins} database route failed.')
            return 0, 0
    except:
        print(f'Invalid species name, cannot be found in {ins} alternate directory.')
        return 0, 0
    
'''
get_annotated_seq function:
    main function that searchs Genbank and Refseq assembly summary files to grab
    relevant accession information for an inputted species, and attempts alternative searches
    if no recent downloads are available from the assembly file
'''
def get_annotated_seq(input_search):
    ##################
    ### Refseq search
    
    input_search = input_search.replace('_', ' ')
    
    print('--------------------------')
    print(f'First searching Refseq for {input_search}...')
    
    if not os.path.exists('absense_fasta_downloads'):
        !mkdir ./absense_fasta_downloads/
    
    # Download the file if it does not exist
    if os.path.exists('./absense_fasta_downloads/assembly_summary_refseq.txt'):
        print('Refseq assembly summary already downloaded! Grabbing species now.')
    else:
        print('Grabbing Refseq assembly from online...')
        # Download file
        try:
            wget.download("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", './absense_fasta_downloads/')
        except:
            print('Error downloading assembly_summary_refseq.txt from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/. Please try manually downloading into the absense_fasta_downloads directory if the problem persists.')
    data_assembly = pd.read_csv('./absense_fasta_downloads/assembly_summary_refseq.txt', sep='\t', skiprows=1)
    
    under_input_search = input_search.replace(' ', '_')
    
    # Grab link
    lnk, acc_name = get_link(data_assembly, input_search)
    
    # Ensure no spaces
    if acc_name != 0:
        acc_name = acc_name.replace(' ', '_')
    
    if lnk != 0 and acc_name != 0:
        # Check if annotation is present
        if os.path.exists('./absense_fasta_downloads/' + acc_name + '_protein.faa'):
            print('Refseq file called ' + acc_name + '_protein.faa already in your directory!')
            return input_search, acc_name
        elif os.path.exists('./absense_fasta_downloads/' + acc_name + '_protein.faa.gz'):
            print('Zipped Refseq file called ' + acc_name + '_protein.faa already in your directory! Unzipping now.')
            zip_dir = acc_name + '_protein.faa.gz'
            !gunzip ./absense_fasta_downloads/$zip_dir
            return input_search, acc_name
        else:
            try:
                wget.download(lnk, './absense_fasta_downloads/')
                zip_dir = acc_name + '_protein.faa.gz'
                !gunzip ./absense_fasta_downloads/$zip_dir
                print('Successfully downloaded and unzipped Refseq file called ' + acc_name + '_protein.faa (for species ' + input_search + ')')
                return input_search, acc_name
            except:
                print('No Refseq annotation exists in the directory for ' + acc_name + ' (' + input_search + ')')
                print('Now searching Genbank...')
    else:
        print('No Refseq annotation exists for ' + input_search)
        print(f'Now searching Genbank for {input_search}...')

    ##################
    ### Genbank search
    
    # Download the file if it does not exist
    if os.path.exists('./absense_fasta_downloads/assembly_summary_genbank.txt'):
        print('Genbank assembly summary already downloaded! Grabbing species now.')
    else:
        print('Grabbing Genbank assembly from online...')
        # Download file
        try:
            wget.download("https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt", './absense_fasta_downloads/')
        except:
            print('Error downloading assembly_summary_refseq.txt from https://ftp.ncbi.nlm.nih.gov/genomes/genbank/. Please try manually downloading into the absense_fasta_downloads directory if the problem persists.')
    data_assembly = pd.read_csv('./absense_fasta_downloads/assembly_summary_genbank.txt', sep='\t', skiprows=1)
    
    # Grab link
    lnk, acc_name = get_link(data_assembly, input_search)
    
    # Ensure no spaces
    if acc_name != 0:
        acc_name = acc_name.replace(' ', '_')
    
    if lnk != 0 and acc_name != 0:
        # Check if annotation is present
        if os.path.exists('./absense_fasta_downloads/' + acc_name + '_protein.faa'):
            print('Genbank file called ' + acc_name + '_protein.faa already in your directory!')
            return input_search, acc_name
        elif os.path.exists('./absense_fasta_downloads/' + acc_name + '_protein.faa.gz'):
            print('Zipped Genbank file called ' + acc_name + '_protein.faa already in your directory! Unzipping now.')
            zip_dir = acc_name + '_protein.faa.gz'
            !gunzip ./absense_fasta_downloads/$zip_dir
            return input_search, acc_name
        else:
            try:
                wget.download(lnk, './absense_fasta_downloads/')
                zip_dir = acc_name + '_protein.faa.gz'
                !gunzip ./absense_fasta_downloads/$zip_dir
                print('Successfully downloaded and unzipped Genbank file called ' + acc_name + '_protein.faa (' + input_search + ')')
                return input_search, acc_name
            except:
                print('No Genbank annotation exists in the directory for ' + acc_name + ' (' + input_search + ')')
    
    ##################
    ### Alternate latest_assembly_version search
    
    # Last ditch effort if there are no species in Refseq or Genbank through the standard way
    # We essentially traverse a different directory and check latest_assembly_versions manually
    print('Trying alternative database route for Refseq and/or Genbank...')
    
    spcies, acc_name = alt_search(input_search, 'REF')
    if acc_name != 0:
        return input_search, acc_name
    
    spcies, acc_name = alt_search(input_search, 'GEN')
    if acc_name != 0:
        return input_search, acc_name
    
    print('Invalid species name or species does not have a protein annotation in either Genbank or Refseq.')
    return 0, 0

'''
get_eukary_odb10 function:
    ensure BUSCO eukaryota genes are locally downloaded for ortholog detection
'''
def get_eukary_odb10():
    # Download the file if it does not exist
    if os.path.exists('eukaryota_odb10.2020-09-10.tar.gz'):
        print('Zipped Eukaryota ODB10 BUSCO genes already downloaded! Unzipping now.')
        !gunzip ./eukaryota_odb10.2020-09-10.tar.gz
        !tar -xvf eukaryota_odb10.2020-09-10.tar
        return True
    elif os.path.exists('./eukaryota_odb10'):
        print('Eukaryota ODB10 BUSCO genes already downloaded and successfully opened.')
        return True
    elif os.path.exists('eukaryota_odb10.2020-09-10.tar'):
        print('Eukaryota ODB10 BUSCO genes already downloaded!')
        !tar -xvf eukaryota_odb10.2020-09-10.tar
        print('Successfully opened tar file.')
        return True
    else:
        print('Grabbing BUSCO genes from its online directory (https://busco-data.ezlab.org/v4/data/lineages/)...')
        
        # Download file using wget
        try:
            print('Downloading...')
            wget.download('https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2020-09-10.tar.gz', './')
            !gunzip ./eukaryota_odb10.2020-09-10.tar.gz
            print('Successfully downloaded.')
            !tar -xvf eukaryota_odb10.2020-09-10.tar
            print('Successfully opened tar file.')
            return True
        except:
            print('Error. File download failed from https://busco-data.ezlab.org/v4/data/lineages/.')
            print('Eukaryotic file name may have changed. Should be eukaryota_odb10.2020-09-10.tar.gz.')
            return False
 
'''
run_hmmer_busco function:
    perform a HMMER homology search with each species protein annotation FASTA file with each
    eukaryota BUSCO gene file and save results in a directory for future analysis
'''
def run_hmmer_busco(species_fasta_list):
    # Create file
    if not os.path.exists('absense_hmmer_out'):
        !mkdir ./absense_hmmer_out/

    # Iterate through species
    for i in species_fasta_list:
        # Go through each .hmm file
        for k in os.listdir('./eukaryota_odb10/hmms'):
            # Ensure proper directory traversal
            dir_name = 'out_hmmer_' + i[0:-4]
            file_name = 'res_hmmer_' + k.strip('.hmm') + '.txt'
            if not os.path.exists(f'./absense_hmmer_out/{dir_name}'):
                !mkdir ./absense_hmmer_out/$dir_name
                
            # Perform hmmsearch
            !hmmsearch ./eukaryota_odb10/hmms/$k ./absense_fasta_downloads/$i > ./absense_hmmer_out/$dir_name/$file_name
            
'''
traverse_hmmer_file function:
    iterate through HMMER functions and grab genes in which match met the inputted E-value threshold
    and also in which there is only 1 mutual hit
'''
def traverse_hmmer_file(file_name, e_val):
    # Iterate through file, check E-val and quantity conditions, and return
    # gene label and query file name
    with open(file_name, "r") as file:
        lines = []
        for line in file:
            lines.append(line)
            
        # Grab E-value and query name
        if (lines[15] == '\n') and (float(lines[14].split()[0]) < e_val):
            return lines[14].split()[8], lines[5].split()[4][23:][0:-4]
        else:
            return 0, 0
   
'''
scan_files function:
    go through HMMER functions and find genes that are consistently hit across all species of interest,
    acting as a heuristic method for ortholog detection with a query of genes of interest to search
'''
def scan_files(species, check_query, e_val, comb):
    # Traverse the relevant files
    dir_name = './absense_hmmer_out/out_hmmer_' + species[0:-4]
    
    gene = []
    query = []
    
    # If check_query is empty, that means it is the cntrl species
    if len(check_query) == 0:
        for f in os.listdir(dir_name):
            g, q = traverse_hmmer_file(dir_name + '/' + f, e_val)
            if g != 0:
                gene.append(g)
            if q != 0:
                query.append(q)
            if g != 0 and q != 0:
                comb[q] = []
                comb[q].append(g)
        return gene, query, comb
    # If check_query is nonempty, then we want to compare only query files that
    # successfully identified an ortholog
    else:
        for f in check_query:
            g, q = traverse_hmmer_file(dir_name + '/' + 'res_hmmer_' + f + '.txt', e_val)
            if g != 0:
                gene.append(g)
            if q != 0:
                query.append(q)
            if g != 0 and q != 0:
                comb[q].append(g)
        return gene, query, comb
    
'''
grab_orthologs function:
    starting with the source species, find all genes that were hit of the eukaryotic genes, and
    repeat the process for all subsequent target species, only search genes that are orthologs of all
    the species searched prior
'''
def grab_orthologs(all_species, e_val):
    # Iterate over species to determine ortholog genes
    
    species_genes = ['']
    query_genes = []
    comb = {}
    
    cntrl_species = all_species[0]
    other_species = all_species[1:]
    
    # Traverse the relevant files
    species_genes[0], query, comb = scan_files(cntrl_species, [], e_val, comb)
    query_genes.append(set(query))
    
    if len(other_species) > 0:
        # Go over other species and widdle down relevant ortholog genes
        for k in other_species:
            tmp, query, comb = scan_files(k, query, e_val, comb)
            
            species_genes.append(tmp)
            query_genes.append(set(query))

    # Take intersection of the set to get only genes that are consistent across species
    out_query = set.intersection(*query_genes)

    return_vals = []
    for i in out_query:
        return_vals.append(comb[i])
        
    # Return an array, where each index is a list of the ortholog genes
    return return_vals

'''
create_muscle_files function:
    iterate through protein annotation files and grab all ortholog genes,
    grouping them into their own FASTA files
'''
def create_muscle_files(adj_list, fasta_files, spcies_list):
    transposed_list = np.array(adj_list).T
    
    # Create directory
    dir_nm = 'absense_files_for_muscle_' + fasta_files[0][0:-4]
    !rm -rf $dir_nm
    if not os.path.exists(dir_nm):
        !mkdir ./$dir_nm

    # Iterate through ortholog genes 
    orthos = [['']*len(transposed_list)]*len(adj_list)
    for i in range(len(transposed_list)):
        with open('./absense_fasta_downloads/' + fasta_files[i], 'r') as file:
            gene_grab = [False]*len(adj_list)
            # Go through and reformat so we can write to proper files
            for j in range(len(adj_list)):
                for line in file:
                    if transposed_list[i][j] in line:
                        gene_grab[j] = True
                    elif '>' in line:
                        gene_grab[j] = False
                    if gene_grab[j] == True and '>' not in line:
                        orthos[j][i] = orthos[j][i] + line.strip('\n').strip('*')

    # Combine the ortholog genes into individual files
    for k in range(len(transposed_list[0])):
        f_name = f'./{dir_nm}/{transposed_list[0][k]}.faa'
        file = open(f_name,'w')
        file.close()
        
        with open(f_name, 'w') as f:
            len_of_line = 70
            for m in range(len(transposed_list)):
                f.write('>' + spcies_list[m][0] + '_' + spcies_list[m].replace(' ', '_').split('_')[1][0:5] + adj_list[k][m] + '\n')

                # Write with proper line breaks
                for p in range(0, len(orthos[k][m]), len_of_line):
                    f.write(orthos[k][m][p:(p+len_of_line)])
                    
                    if p == len(orthos[k][m]) - (len(orthos[k][m]) % 70):
                        f.write('\n')

'''
compute_muscle function:
    going through the created FASTA files for MUSCLE, perform MUSCLE locally on each
    set of ortholog genes and store the multiple sequence alignment results
'''
def compute_muscle(accessions):
    source_accession = accessions[0]
    
    # Directory in which the data resides
    dirloc = f"./absense_files_for_muscle_{source_accession[0:-4]}" 
    fl_str = f"absense_files_for_muscle_{source_accession[0:-4]}"

    count = 0
    lst_file = []
    # Check if valid files
    for file in os.scandir(dirloc):
        if (file.path.endswith(".faa")) and file.is_file():
            count = count + 1
            lst_file.append(file.name)
     
    # Ensure correct directory exists
    fl_name = './muscle_out_for_absense_' + fl_str[25:]
    if not os.path.exists(fl_name):
        !mkdir ./$fl_name

    # Assumes the current directory the script resides in contains the 'muscle3.8.31_i86win32.exe' file
    for i in range(count):
        # Local version
        in_dir = fl_str
        out_dir = fl_name
        val = lst_file[i]
        out_nm = 'out_' + lst_file[i][0:9] + '.faa'
        in_file = f'{in_dir}/{val}'
        out_file = f'{out_dir}/{out_nm}.faa'

        !muscle -in ./$in_dir/$val -out ./$out_dir/$out_nm
        
        # API Version (for possible future reference)
        # os.system('cmd /k "python muscle.py --email {email} {val}"'.format(val = lst_file[i-1]))
    
'''
compute_protdist function:
    from the multiple sequence alignment files, run Protdist from PHYLIP
    on the MSA files to get the evolutionary distances
'''
def compute_protdist(source_accession):
    # Directory in which the data resides
    dirloc = f"./muscle_out_for_absense_{source_accession[0:-4]}" 
    fl_str = f'muscle_out_for_absense_{source_accession[0:-4]}'

    lst_file = []

    for file in os.scandir(dirloc):
        if (file.path.endswith(".faa")) and file.is_file():
            lst_file.append(file)

    # Init dictionary for gene names and respective sequences
    genes = {}

    # Read first file to get the gene names
    with open(lst_file[0], 'r') as t_file:
        for l in t_file:
            if l[0] == '>':
                genes[l[1:8]] = ''

    # Iterate through the file contents
    for i in range(len(lst_file)):
        # Open file to read it
        gene_track = ''
        with open(lst_file[i], 'r') as r_file:
            for line in r_file:
                if line[0] == '>':
                    gene_track = line[1:8]
                else:
                    genes[gene_track] = genes[gene_track] + line.strip('\n').strip('*')

    fl_name = './protdist_infile_outfile_' + fl_str[23:]
    if not os.path.exists(fl_name):
        !mkdir ./$fl_name

    if os.path.exists('outfile'):
        !rm outfile

    # Write and format file for PHYLIP and Protdist
    # Must have number of species and length of each sequence at the top left, followed by each name followed by a sequence
    # *Has* to be named infile (as per Protdist documentation, seems like there's no way to change it)
    with open('./infile', 'w') as w_file:
        w_file.write('   ' + str(len(genes.keys())) + '   ' + str(len(genes[list(genes.keys())[0]])) + '\n')
        for k in genes:
            w_file.write(k + '   ' + genes[k] + '\n')

    !echo Y | ./protdist

    !cp infile infile ./$fl_name
    !cp outfile outfile ./$fl_name
    
    return fl_name

'''
read_accession_file function:
    read the accession file of all the species of interest and their corresponding 
    protein annotation FASTA file names for use in the pipeline
'''
def read_accession_file(file_name):
    # Go through accession file and pass back the accession codes and species names
    species_names = []
    species_accessions = []
    
    # Iterate file and return valid accessions and species names
    with open(file_name, 'r') as f:
        k = 0
        for line in f:
            if k != 0:
                # Split comma delineated file
                txt = line.strip('\n').strip(' ').split(',')
                if txt[1] != '0_protein.faa':
                    species_names.append(txt[0])
                    species_accessions.append(txt[1])
            k += 1
    
    return species_names, species_accessions

'''
grab_assemblies function:
    iterate through species of interest and perform an online depth search through Refseq and Genbank
    to download the protein annotations to a local directory
'''
def grab_assemblies(list_of_species, override_dict = {}):
    # List of species contains list of species, where first species is the source and the others are the target
    # override_dict is a dictionary of species and their corresponding FASTA files, for when you manually input
    # FASTA files into the directory so the pipeline can incorporate them correctly
    
    species_names = [0]*len(list_of_species)
    species_accessions = [0]*len(list_of_species)
    
    # Grab protein annotations for all species
    for i in range(len(list_of_species)):
        species_names[i], species_accessions[i] = get_annotated_seq(list_of_species[i])
        
    # Write output accession codes to a text file
    with open('absense_accession_codes_infile.txt', 'w') as f:
        f.write('name,accession,type\n')
        for k in range(len(species_names)):
            if list_of_species[k] in override_dict:
                f.write(f'{list_of_species[k]},{override_dict[list_of_species[k]][0]},{override_dict[list_of_species[k]][1]}\n')
            elif k != 0:
                f.write(f'{list_of_species[k]},{species_accessions[k]}_protein.faa,target\n')
            else:
                f.write(f'{list_of_species[k]},{species_accessions[k]}_protein.faa,source\n')
            
    f.close()
    
    if 0 in species_accessions:
        print('--------------------')
        print('Warning: not all species had their protein annotations downloaded.')
        print('Try manually downloading the necessary files and adding the species names and protein annotation file names to')
        print('the absense_accession_codes_infile.txt file (or something else if it was changed).')
        print('When satisfied, proceed with the evolutionary distance computation.')
    else:
        print('--------------------')
        print('Successfully downloaded protein annotations for all species!')
     
'''
run_evo_dist_computation function:
    perform the pipeline of getting the eukaryotic BUSCO genes, performing HMMER homology search
    to detect orthologs, getting the MSA of the orthologous genes, and then calculating the distances
    based on the alignment
'''
def run_evo_dist_computation(evalue = 0.01, acc_file = 'absense_accession_codes_infile.txt'):
    # Ensure e-value is below HMMER threshold
    if evalue > 0.01:
        print('Warning. E-value is above HMMER threshold of 0.01. Please reduce for meaningful results.')
        return
    
    # Read accession codes file
    print('-------------------')
    species_names, species_accessions = read_accession_file(acc_file)
        
    # Grab BUSCO genes
    success_busco = get_eukary_odb10()
    
    if success_busco == False:
        print('Error downloading BUSCO genes from https://busco-data.ezlab.org/v4/data/lineages/')
        print('Please try manually downloading into the local directory if the issue persists.')
        print('Program cannot continue without the BUSCO eukaryota odb10 genes.')
        return False
    
    print('-------------------')
    print('Running HMMER calcuations now...')
    run_hmmer_busco(species_accessions)
    
    print('-------------------')
    print(f'Iterating results and grabbing orthologs now with E-value of {evalue}...')
    result = grab_orthologs(species_accessions, evalue)
    
    print('-------------------')
    print('Creating files for MUSCLE now...')
    create_muscle_files(result, species_accessions, species_names)
    
    print('-------------------')
    print('Running MUSCLE now...')
    compute_muscle(species_accessions)
    
    print('-------------------')
    print('Running Protdist now...')
    out = compute_protdist(species_accessions[0])
    
    print('-------------------')
    print(f'Done! Results are in "outfile" in the {out} directory. Thanks!')

### cell
grab_assemblies(list_of_species = 
                ['Saccharomyces cerevisiae',
                 'Saccharomyces mikatae',
                 'Naumovozyma castellii',
                 'Kluyveromyces lactis',
                 'Aspergillus nidulans',
                 'Saccharomyces kudriavzevii',
                 'Yarrowia lipolytica',
                 'Schizosaccharomyces pombe',
                 'Eremothecium gossypii',
                 'Saccharomyces paradoxus',
                 'Saccharomyces bayanus'],
                 override_dict = 
                 {'Saccharomyces mikatae' : ['Smik.faa' , 'target'],
                  'Saccharomyces bayanus' : ['Sbay.faa', 'target']})

### cell
run_evo_dist_computation(evalue = 0.01, acc_file = 'absense_accession_codes_infile.txt')
