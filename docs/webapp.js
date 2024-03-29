const DEFAULT_MSG = "Interpret undetected homologs. Runs entirely in your browser using React, MaterialUI, and Pyodide. Thanks to scikit-hep/repo-review for inspiration!";

function Heading(props) {
  return (
    <MaterialUI.Box sx={{ flexGrow: 1, mb: 2 }}>
      <MaterialUI.AppBar position="static">
        <MaterialUI.Toolbar>
          <MaterialUI.Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
            Try abSENSE today!
          </MaterialUI.Typography>
          <MaterialUI.Button
            href="https://github.com/PrincetonUniversity/abSENSE"
            color="inherit">Source</MaterialUI.Button>
          <MaterialUI.Button
            href="https://github.com/scikit-hep/repo-review"
            color="inherit">Inspiration</MaterialUI.Button>
        </MaterialUI.Toolbar>
      </MaterialUI.AppBar>
    </MaterialUI.Box>
  );
}

function Results(props) {
  const output = [];

  if(typeof props.results === 'string'){
    return (
      <MaterialUI.Box sx={{bgcolor: 'background.paper', width: 1}}>
        <img
          src={`${props.results}`}
          srcSet={`${props.results} 1x`}
          alt="Plot of the gene"
          loading="lazy"
          style={{maxWidth: "100%"}}
        />
      </MaterialUI.Box>
    )
  }

  props.results.forEach((value, name) => {
    const results_components = value.map(result => {
      const icon = <MaterialUI.Icon>{name == 'channels' ? 'looks_one' : (name == 'samples' ? 'looks_two' : 'looks_3')}</MaterialUI.Icon>;
      const msg = (
        <React.Fragment>
          <MaterialUI.Typography
            sx={{ display: 'inline' }}
            component="span"
            variant="body2"
            color="text.primary"
          >
          {JSON.stringify(result)}
          </MaterialUI.Typography>
        </React.Fragment>
      )


      return (
        <MaterialUI.ListItem disablePadding key='{name}-{value}'>
          <MaterialUI.ListItemButton>
            <MaterialUI.ListItemIcon>
              {icon}
            </MaterialUI.ListItemIcon>
            <MaterialUI.ListItemText primary={msg} />
          </MaterialUI.ListItemButton>
        </MaterialUI.ListItem>
      );
    });

    output.push(
      <li key={`section-${name}`}>
        <ul>
          <MaterialUI.ListSubheader>{name}</MaterialUI.ListSubheader>
          {results_components}
        </ul>
      </li>
    );
  });

  return (
    <MaterialUI.Box sx={{bgcolor: 'background.paper'}} >
      <MaterialUI.List subheader={<li />} overflow='auto'>{output}</MaterialUI.List>
    </MaterialUI.Box>
  );
}

async function prepare_pyodide(app) {
  const pyodide = await loadPyodide({
    stdout: (message) => {app.setState({msg: message, progress: true})}
  });

  await pyodide.loadPackage('micropip', (message) => {app.setState({msg: message, progress: true})});
  app.setState({msg: "Loading abSENSE packages ... ", progress: true});
  await pyodide.runPythonAsync(`
    import micropip
    await micropip.install(["setuptools", "matplotlib", "absense==0.1.6"])
    import io
    from abSENSE import pyodide_utilities
    analyzer = None
    missing_homologs = []
  `);
  app.setState({msg: DEFAULT_MSG, progress: false});
  app.setState({absense_ready: true})
  return pyodide;
}

function MyThemeProvider(props) {
  const prefersDarkMode = MaterialUI.useMediaQuery('(prefers-color-scheme: dark)');

  const theme = React.useMemo(
    () =>
    MaterialUI.createTheme({
      palette: {
        mode: prefersDarkMode ? 'dark' : 'light',
      },
    }),
    [prefersDarkMode],
  );

  return (
    <MaterialUI.ThemeProvider theme={theme}>
    {props.children}
    </MaterialUI.ThemeProvider>
  );
}

function a11yPropsTab(index) {
  return {
    id: `simple-tab-${index}`,
    'aria-controls': `simple-tabpanel-${index}`,
  };
}

function TabPanel(props) {
  const { children, value, index, ...other } = props;

  return (
    <div
    role="tabpanel"
    hidden={value !== index}
    id={`simple-tabpanel-${index}`}
    aria-labelledby={`simple-tab-${index}`}
    {...other}
    >
    {value === index && (
      <MaterialUI.Box sx={{ p: 3 }}>
      {children}
      </MaterialUI.Box>
    )}
    </div>
  );
}

TabPanel.propTypes = {
  children: PropTypes.node,
  index: PropTypes.number.isRequired,
  value: PropTypes.number.isRequired,
};

class App extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      results: new Map(),
      text_input: `Species	Distance	NP_010181.2	NP_009362.1	NP_014555.1	NP_116682.3	NP_011284.1	NP_011878.1	NP_013320.1	NP_014160.2	NP_014890.1
S_cer	0.000	2284	1027	712	352	308	588	1491	941	394
S_par	0.051	2254	1016	698	171	225	522	1205	927	264
S_mik	0.088	2234	1018	696	92.8	179	446	1013	889	254
S_kud	0.104	2197	1010	688	82	178	455	996	893	271
S_bay	0.108	2200	1008	684	0	191	439	910	885	230
S_cas	0.363	1895	909	N/A	0	0	0	141	638	58.9
K_wal	0.494	1767	903	617	0	0	0	101	0	0
A_gos	0.518	1736	875	N/A	0	0	0	63.2	657	0
K_lac	0.557	1691	894	611	0	0	0	102	0	0
A_nid	0.903	1099	691	N/A	0	0	0	0	378	0
S_pom	0.922	1194	684	344	0	0	0	0	0	0
Y_lip	0.954	1297	714	N/A	0	0	0	0	375	0`,
      distanceString: null,
      distanceFile: null,
      bitscoreString: null,
      bitscoreFile: null,
      msg: "Initializing ... ",
      progress: false,
      tab: "textInput",
      e_value_cutoff: 0.001,
      gene_length: 400,
      db_size: 4000000,
      gene_select_name: null,
      gene_names: [],
      missing_homologs: [],
      absense_ready: false
    };
  }
  componentDidMount() {
    this.pyodide_promise = prepare_pyodide(this);
  }

  handleFileChangeDistance = (e) => {
    const distanceFile = e.target.files[0];
    if (!distanceFile) return;

    const reader = new FileReader();
    reader.onloadend = e => {
      this.setState({distanceString: e.target.result, progress: false, msg: `Loaded ${this.state.distanceFile.name}`});
    }
    reader.readAsText(distanceFile);
    this.setState({distanceFile: distanceFile});
  }

  handleFileChangeBitscore = (e) => {
    const bitscoreFile = e.target.files[0];
    if (!bitscoreFile) return;

    const reader = new FileReader();
    reader.onloadend = e => {
      this.setState({bitscoreString: e.target.result, progress: false, msg: `Loaded ${this.state.bitscoreFile.name}`});
    }
    reader.readAsText(bitscoreFile);
    this.setState({bitscoreFile: bitscoreFile});
  }

  handleChangeMultiple = (e) => {
    this.state.gene_select_name = e.target.value;
    this.updatePlot();
  }

  updatePlot(){
    if(this.state.gene_select_name && this.state.gene_names){
      this.pyodide_promise.then(pyodide => {
        try{
          let results = pyodide.runPython(`
                                print(analyzer.bitscores.loc["${this.state.gene_select_name}"])
                                plot, missing_homologs = pyodide_utilities.get_plot_and_table("${this.state.gene_select_name}", analyzer)
                                img_str = pyodide_utilities.generate_png_string(plot)
                                img_str
                            `);
          this.setState({
            results: results,
            missing_homologs: pyodide.globals.get('missing_homologs').toJs(),
            msg: `Plot for ${this.state.gene_select_name}`,
            progress: false
          });
        }
        catch(e){
          console.log(e.message)
          this.setState({msg: `Error during analysis, full traceback in debug console`,
                         progress: false});
        }
      });
    }
  }

  async handlePlotUploadFile() {
    if (!this.state.bitscoreString && !this.state.distanceString) {
      alert(`Please provide values in distance and bitscore files!`);
      return;
    }
    this.setState({"results": new Map(), msg: "Running...", progress: true});
    this.pyodide_promise.then(async pyodide => {
      try {
        window.bitscore = this.state.bitscoreString;
        window.distance = this.state.distanceString;
        let gene_list = await pyodide.runPythonAsync(`
            from js import bitscore, distance
            bitscore = io.StringIO(bitscore)
            distance = io.StringIO(distance)
            e_value_cutoff = float(${this.state.e_value_cutoff})
            gene_length = float(${this.state.gene_length})
            db_size = float(${this.state.db_size})
            analyzer = pyodide_utilities.get_analyzer_from_files(bitscore, distance, e_value_cutoff, gene_length, db_size)
            print('Done analyzing')
            list(analyzer.genes)
        `);
        gene_list = gene_list.toJs();
        delete window.bitscore;
        delete window.distance;
        if(gene_list){
          this.state.gene_names = gene_list
          this.state.gene_select_name = gene_list[0]
          this.updatePlot();
        }
        else
          this.setState({msg: `Error during analysis`,
            progress: false});
      }
      catch(e){
        console.log(e.message)
        this.setState({msg: `Error during analysis, full traceback in debug console`,
          progress: false});
      }
    });
  }

  async handlePlotTextInput() {
    if (!this.state.text_input) {
      alert(`Please provide score table!`);
      return;
    }
    this.setState({"results": new Map(), msg: "Running...", progress: true});

    this.pyodide_promise.then(async pyodide => {
      try {
        window.text_input = this.state.text_input;
        let gene_list = await pyodide.runPythonAsync(`
                                from js import text_input
                                e_value_cutoff = float(${this.state.e_value_cutoff})
                                gene_length = float(${this.state.gene_length})
                                db_size = float(${this.state.db_size})
                                analyzer = pyodide_utilities.get_analyzer_from_text(text_input, e_value_cutoff, gene_length, db_size)
                                print('Done analyzing')
                                list(analyzer.genes)
                            `);
        gene_list = gene_list.toJs();
        delete window.text_input;
        if(gene_list){
          this.state.gene_names = gene_list
          this.state.gene_select_name = gene_list[0]
          this.updatePlot();
        }
        else
          this.setState({msg: `Error during analysis`,
            progress: false});
      }
      catch(e){
        console.log(e.message)
        this.setState({msg: `Error during analysis, full traceback in debug console`,
          progress: false});
      }
    });
  }


  async handleAnalyze() {
    if (this.state.tab == "textInput")
      return this.handlePlotTextInput();

    if (this.state.tab == "publishedFungi"){
      let gitdata = await fetch("https://raw.githubusercontent.com/caraweisman/abSENSE/c355c458e83722a0ffdf7284d4ea1f6f29ce205f/Fungi_Data/Fungi_Bitscores");
      this.state.bitscoreString = await gitdata.text();
      gitdata = await fetch("https://raw.githubusercontent.com/caraweisman/abSENSE/c355c458e83722a0ffdf7284d4ea1f6f29ce205f/Fungi_Data/Fungi_Distances");
      this.state.distanceString = await gitdata.text();
    }

    if (this.state.tab == "publishedInsect"){
      let gitdata = await fetch("https://raw.githubusercontent.com/caraweisman/abSENSE/c355c458e83722a0ffdf7284d4ea1f6f29ce205f/Insect_Data/Insect_Bitscores");
      this.state.bitscoreString = await gitdata.text();
      gitdata = await fetch("https://raw.githubusercontent.com/caraweisman/abSENSE/c355c458e83722a0ffdf7284d4ea1f6f29ce205f/Insect_Data/Insect_Distances");
      this.state.distanceString = await gitdata.text();
    }

    return this.handlePlotUploadFile();
  }

  savePlot() {
    if(this.state.gene_select_name && this.state.gene_names){
      this.pyodide_promise.then(pyodide => {
        let output = pyodide.runPython(`
                            from js import Blob
                            plot, missing_homologs = pyodide_utilities.get_plot_and_table("${this.state.gene_select_name}", analyzer)
                            buf = io.BytesIO()
                            plot.save(buf, format='svg')
                            buf.seek(0)
                            buf.read().decode('UTF-8')
                        `);
        var blob = new Blob([output], {type: 'application/text'});
        var elem = window.document.createElement('a');
        elem.href = window.URL.createObjectURL(blob)    ;
        elem.download = `${this.state.gene_select_name}.svg`;
        document.body.appendChild(elem);
        elem.click();        
        document.body.removeChild(elem);

      });
    }
  }

  render() {

    let icon_close_distance = ''
    if(this.state.distanceFile){
      icon_close_distance = (
        <MaterialUI.IconButton
        aria-label="close"
        onClick={() => {this.setState({distanceFile: null, text_input: null, msg: `Deleted ${this.state.distanceFile.name}`})}}
        sx={{
          color: (theme) => theme.palette.grey[500],
        }}
        >
        <MaterialUI.Icon>close</MaterialUI.Icon>
        </MaterialUI.IconButton>
      )
    }

    let icon_close_bitscore = ''
    if(this.state.bitscoreFile){
      icon_close_bitscore = (
        <MaterialUI.IconButton
        aria-label="close"
        onClick={() => {this.setState({bitscoreFile: null, text_input: null, msg: `Deleted ${this.state.bitscoreFile.name}`})}}
        sx={{
          color: (theme) => theme.palette.grey[500],
        }}
        >
        <MaterialUI.Icon>close</MaterialUI.Icon>
        </MaterialUI.IconButton>
      )
    }
    let select_gene_window = ''
    if(this.state.gene_select_name){
      select_gene_window = (
        <MaterialUI.Stack alignItems="center" sx={{ m: 1, mb: 0 }}>
          <MaterialUI.Box>
            <MaterialUI.Typography variant="button" display="block" gutterBottom>
              Select a Gene
            </MaterialUI.Typography>
          </MaterialUI.Box>
          <MaterialUI.Typography variant="button" display="block" gutterBottom>
            <MaterialUI.Select
            multiple
            native
            variant="filled"
            value={this.gene_select_name}
            onChange={this.handleChangeMultiple}
            label="List of Genes"
            inputProps={{size: 30}}>
              {this.state.gene_names.map((name, i) => (i==0) ?
                (<option key={name} value={name} selected>{name}</option>) :
                (<option key={name} value={name}>{name}</option>)
              )}
            </MaterialUI.Select>
          </MaterialUI.Typography>
        </MaterialUI.Stack>
      )
    }

    return !this.state.absense_ready ? (
      <div className="absense_load_progress">
        <MaterialUI.Typography variant="h2">Loading abSENSE</MaterialUI.Typography>
        <MaterialUI.Typography variant="body1" component="div">
          {this.state.msg}
        </MaterialUI.Typography>
        <MaterialUI.LinearProgress disabled={!this.state.progress} />
      </div>
    ) : (
      <MyThemeProvider>
        <MaterialUI.CssBaseline />
        <MaterialUI.Box>
          <Heading />
          <MaterialUI.Grid container>
            <MaterialUI.Grid item xs={2}>
              <MaterialUI.Stack spacing={2} alignItems="left" sx={{ m: 0, mb: 0 }} >
                <MaterialUI.Tabs sx={{flexGrow: 2}} orientation='vertical' value={this.state.tab} onChange={(e, v) => {this.setState({tab: v})}} aria-label="basic tabs example">
                  <MaterialUI.Tab label="Text Input" value="textInput" />
                  <MaterialUI.Tab label="File Upload" value="uploadFile" />
                  <MaterialUI.Tab label="Fungi Results" value="publishedFungi" />
                  <MaterialUI.Tab label="Insect Results" value="publishedInsect" />
                </MaterialUI.Tabs>
                <MaterialUI.Button onClick={() => this.handleAnalyze()}
                  variant="contained"
                  size="large"
                  disabled={this.state.progress || (this.state.tab == "uploadFile" && !(this.state.bitscoreFile && this.state.distanceFile)) || (!this.state.text_input && this.state.tab == "textInput")}>Analyze</MaterialUI.Button>
                <MaterialUI.Button onClick={() => this.savePlot()} variant="contained" size="large" disabled={!(this.state.gene_select_name && this.state.gene_names)}>Save Plot</MaterialUI.Button>
              </MaterialUI.Stack>
            </MaterialUI.Grid>
            <MaterialUI.Grid item xs={8}>
              <TabPanel value={this.state.tab} index="textInput">
                <MaterialUI.TextField 
                  id="workspace-text"
                  label="Custom Analysis"
                  helperText="Separate columns with pipe(|), comma(,) or tab."
                  variant="filled"
                  autoFocus={true}
                  sx={{flexGrow: 1, m: 1}}
                  value={this.state.text_input}
                  multiline
                  fullWidth
                  minRows={3}
                  maxRows={9}
                  onInput={(e) => this.setState({text_input: e.target.value})}
                  />
              </TabPanel>
              <TabPanel value={this.state.tab} index="uploadFile">
                <MaterialUI.Stack alignItems="left" sx={{ m: 0, mb: 0 }} spacing={2}>
                  <MaterialUI.Stack direction="row" spacing={2} alignItems="center" sx={{ m: 0, mb: 0 }} >
                    <MaterialUI.Button variant="contained" component="label" disabled={this.state.progress}>
                      Upload Distance <input type="file" onChange={this.handleFileChangeDistance} value={!this.state.distanceFile ? '' : null} hidden />
                    </MaterialUI.Button>
                    {icon_close_distance}
                    <MaterialUI.Typography variant="body1" component="div" sx={{flexGrow: 2}}>
                      {this.state.distanceFile ? this.state.distanceFile.name : 'Select a distance file.'}
                    </MaterialUI.Typography>
                  </MaterialUI.Stack>
                  <MaterialUI.Stack direction="row" spacing={2} alignItems="center" sx={{ m: 0, mb: 0 }} >
                    <MaterialUI.Button variant="contained" component="label" disabled={this.state.progress}>
                      Upload Bitscore <input type="file" onChange={this.handleFileChangeBitscore} value={!this.state.bitscoreFile ? '' : null} hidden />
                    </MaterialUI.Button>
                    {icon_close_bitscore}
                    <MaterialUI.Typography variant="body1" component="div" sx={{flexGrow: 2}}>
                      {this.state.bitscoreFile ? this.state.bitscoreFile.name : 'Select a bitscore file.'}
                    </MaterialUI.Typography>
                  </MaterialUI.Stack>
                </MaterialUI.Stack>
              </TabPanel>
              <TabPanel value={this.state.tab} index="publishedFungi">
                From <a href="https://github.com/caraweisman/abSENSE/tree/c355c458e83722a0ffdf7284d4ea1f6f29ce205f/Fungi_Data">Fungi_Data</a> as described in <a href="https://doi.org/10.1371/journal.pbio.3000862">Weisman et. al.</a>
              </TabPanel>
              <TabPanel value={this.state.tab} index="publishedInsect">
                From <a href="https://github.com/caraweisman/abSENSE/tree/c355c458e83722a0ffdf7284d4ea1f6f29ce205f/Insect_Data">Insect_Data</a> as described in <a href="https://doi.org/10.1371/journal.pbio.3000862">Weisman et. al.</a>
              </TabPanel>
            </MaterialUI.Grid>
            <MaterialUI.Grid item xs={2}>
              <MaterialUI.Stack alignItems="left" sx={{ m: 0, mb: 0 }} >
                <MaterialUI.TextField
                id="e-value-cutoff"
                label="E-value cutoff"
                variant="filled"
                autoFocus={false}
                sx={{flexGrow: 1, m: 1, width:{sm:200,md:300}}}
                value={this.state.e_value_cutoff}
                onInput={(e) => this.setState({e_value_cutoff: e.target.value})}
                type="number"
                />
                <MaterialUI.TextField
                id="gene-length"
                label="Gene length (aa)"
                variant="filled"
                autoFocus={false}
                sx={{flexGrow: 1, m: 1, width:{sm:200,md:300}}}
                value={this.state.gene_length}
                onInput={(e) => this.setState({gene_length: e.target.value})}
                type="number"
                />
                <MaterialUI.TextField
                id="db-size"
                label="Database size (per species) (aa)"
                variant="filled"
                autoFocus={false}
                sx={{flexGrow: 1, m: 1, width:{sm:200,md:300}}}
                value={this.state.db_size}
                onInput={(e) => this.setState({db_size: e.target.value})}
                type="number"
                />
              </MaterialUI.Stack>
            </MaterialUI.Grid>
            <MaterialUI.Grid item xs={12}>
              <MaterialUI.Stack direction="row" spacing={2} alignItems="left" sx={{ m: 0, mb: 0 }} >
                <MaterialUI.Grid item xs={2}>
                  {select_gene_window}
                </MaterialUI.Grid>
                <MaterialUI.Grid item xs={6}>
                  <MaterialUI.Paper elevation={3}>
                    <Results results={this.state.results} />
                    <MaterialUI.Box sx={{ p: 2 }}>
                      {this.state.progress && <MaterialUI.LinearProgress />}
                      <MaterialUI.Typography variant="body1" component="div">
                      {this.state.msg}
                      </MaterialUI.Typography>
                    </MaterialUI.Box>
                  </MaterialUI.Paper>
                </MaterialUI.Grid>
                <MaterialUI.Grid item xs={4}>
                  <MaterialUI.TableContainer component={MaterialUI.Paper}>
                    <MaterialUI.Table stickyHeader aria-label="Missing homologs" size="small" >
                      <MaterialUI.TableHead>
                        <MaterialUI.TableRow>
                          <MaterialUI.TableCell size="small">
                            Species
                          </MaterialUI.TableCell>
                          <MaterialUI.TableCell  size="small">
                            P(detection<br/>failure | E={this.state.e_value_cutoff})
                          </MaterialUI.TableCell>
                          <MaterialUI.TableCell  size="small">
                            Predicted<br/>bitscore
                          </MaterialUI.TableCell>
                          <MaterialUI.TableCell  size="small">
                            99% CI
                          </MaterialUI.TableCell>
                        </MaterialUI.TableRow>
                      </MaterialUI.TableHead>
                      <MaterialUI.TableBody>
                        {this.state.missing_homologs.map((row) => (
                          <MaterialUI.TableRow key={row.get('species')}
                          sx = {{ '&:last-child td, &:last-child th': { border: 0 } }}
                          >
                            <MaterialUI.TableCell component="th" scope="row"  size="small">{row.get('species')}
                            </MaterialUI.TableCell>
                            <MaterialUI.TableCell  size="small">{row.get('p_value')}</MaterialUI.TableCell>
                            <MaterialUI.TableCell  size="small">{row.get('bitscore')}</MaterialUI.TableCell>
                            <MaterialUI.TableCell  size="small">{row.get('interval')}</MaterialUI.TableCell>
                          </MaterialUI.TableRow>
                        ))}
                      </MaterialUI.TableBody>
                    </MaterialUI.Table>
                  </MaterialUI.TableContainer>
                </MaterialUI.Grid>
              </MaterialUI.Stack>
            </MaterialUI.Grid>
          </MaterialUI.Grid>
        </MaterialUI.Box>
      </MyThemeProvider>
    );
  }
}
