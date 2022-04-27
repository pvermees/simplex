// 1. Startup

var glob = {
    "simplex": {
	"method": { "method": ["IGG-UPb"] },
	"samples": null,
	"outliers": null,
	"calibration": null,
	"calibrated": null
    },
    "i": null,
    "start": true,
    "names": null,
    "class": ["simplex"],
    "multi": null,
    "selected": 0,
    "ratios": true,
    "log": true,
    "cov": false,
    "exterr": false,
    "xy": false,
    "shownum": false,
    "calibration": {
	"caltype": null,
	"standcomp": null,
	"standtype": null,
	"preset": null,
	"tst": null,
	"del": {
	    "ratios": null,
	    "delval": null,
	    "delcov": null,
	    "refval": null,
	    "refcov": null
	}
    },
    "sampleprefix": [''],
    "delta": {
	"type":"delta-prime",
	"preset": null,
	"ratios": null,
	"val": null
    },
    "IsoplotRtype":"U-Pb",
    "standards": [],
    "samples": [],
    "buttonIDs": ["setup","drift","logratios","calibration","samples","finish"]
}

function start() {
    // This will be moved into shinylight.initialize()
    new Promise((resolve, reject) => {
        rrpc.initialize(
	    () => resolve(),
	    error => reject(error)
	);
    }).then(() =>
        welcome()
    );
}

function selectedButton() {
    return document.getElementById(
        glob.buttonIDs[glob.selected]
    );
}

function selectButton(i){
    selectedButton().classList.remove('on')
    glob.selected = i;
    selectedButton().classList.add('on')
}

async function loadPage(url){
    let response = await fetch(url);
    let text = await response.text();
    document.getElementById("contents").innerHTML = text;
}

function togglehelp(on,plot=true){
    document.getElementById('help').style.display = on ? 'inline' : 'none';
    if (plot){
	document.getElementById('output').style.display = on ? 'none' : 'inline';
    }
}

async function loadHelp(url,plot=true){
    let response = await fetch(url);
    let text = await response.text();
    togglehelp(true,plot=plot);
    document.getElementById("help").innerHTML = text;
}

function welcome(){
    loadPage("welcome.html");
}

// 2. Setup

async function setup(){
    selectButton(0);
    loadPage("setup.html").then(
        () => {
	    initpreset();
	    if (glob.start) loadPresets()
	    else showPresets()
	    glob.start = false;
	},
        error => alert(error)
    ).then(
	() => loadHelp("help/setup.html",false),
	error => alert(error)
    );
}

function initpreset(){
    document.getElementById("methods").value = glob.simplex.method.method[0];
}

function loadPresets(){
    resetglob();
    const m = document.getElementById("methods").value;
    glob.simplex.method.method[0] = m;
    shinylight.call('presets', { method: m }, null).then(
	result => {
	    result2simplex(result)
	    showPresets();
	    fileFormats();
	},
	err => alert(err)
    )
}

function resetglob(){
    glob.i =  0;
    glob.multi = false; // multi-element
    glob.sampleprefix = [''];
    glob.standards = [];
    glob.samples = [];
    glob.calibration = {
	'caltype': null, // 'average' or 'regression'
	'standcomp': 'preset2stand', // 'manualstand', 'preset2stand', 't2stand' or 'd2stand'
	'standtype': 'measured', // 'measured' or 'commonradio'
	'preset': null, // 'Plesovice', 'NBS28', ...
	'tst': [0,0], 
	'del': {
	    'ratios': null,
	    'val': null,
	    'cov': null,
	    'refval': null,
	    'refcov': null
	},
	'delta': {
	    'type':'delta-prime',
	    'preset': null, // 'VSMOW', 'troilite', ...
	    'ratios': null,
	    'val': null
	},
	'IsoplotRtype':'U-Pb',
    }
}

function method(el){
    glob.simplex.method[el.id] = el.value.split(',');
    glob.class = ['simplex']; // reset calculations
    checkmethod();
}

function renameIons(){
    let ions = glob.simplex.method.ions;
    let keys = Object.keys(glob.simplex.samples);
    for (let i=0; i<keys.length; i++){
	glob.names.samples[keys[i]].detector = ions;
	glob.names.samples[keys[i]].dtype = ions;
	glob.names.samples[keys[i]].dwelltime = ions;
	glob.names.samples[keys[i]].sbm.cnames = ions;
	glob.names.samples[keys[i]].signal.cnames = ions;
	glob.names.samples[keys[i]].time.cnames = ions;
    }
}

function showPresets(){
    let assign = (id) => document.getElementById(id).value =
	glob.simplex.method[id];
    assign('description');
    assign('instrument');
    assign('ions');
    assign('num');
    assign('den');
    assign('bkg');
}

function show(cls){
    let set = document.querySelectorAll(cls);
    if (set==null) return    
    document.querySelectorAll(cls).forEach(
	function(item){
	    item.classList.remove("hidden");
	});
}
function hide(cls){
    let set = document.querySelectorAll(cls);
    if (set==null) return
    document.querySelectorAll(cls).forEach(
	function(item){
	    item.classList.add("hidden");
	});
}

function fileFormats(){
    let accept = ['.asc','.op','.pd'];
    if (glob.simplex.method.instrument=='Cameca'){
	accept = '.asc';
    } else if (glob.simplex.method.instrument=='SHRIMP'){
	accept = ['.op','.pd'];
    } else {
	alert('Unrecognised instrument.')
    }
    document.getElementById('upload').accept = accept;
}

// From https://masteringjs.io/tutorials/fundamentals/filereader
function readFile(file) {
    return new Promise((resolve, reject) => {
	const reader = new FileReader();
	reader.onload = res => {
	    resolve(res.target.result);
	};
	reader.onerror = err => reject(err);
	reader.readAsText(file);
    });
}

// read all files for conversion to textConnection
async function readFiles(){
    let f = document.getElementById('upload').files;
    let fns = {};
    let tcs = {};
    for (let i=0; i<f.length; i++){
	fns[i] = f[i].name;
	tcs[i] = await readFile(f[i]);
    }
    return({fns:fns, tcs:tcs})
}

async function upload(){
    let status = document.getElementById('upload-status');
    readFiles().then(
	f => {
	    status.innerHTML = "<span class='blink'>Reading...</span>";
	    m = glob.simplex.method;
	    shinylight.call('upload', {f:f, m:m}, null).then(
		result => {
		    result2simplex(result);
		    document.getElementById('ions').value =
			glob.simplex.method.ions.toString();
		    document.getElementById('num').value =
			glob.simplex.method.num.toString();
		    document.getElementById('den').value =
			glob.simplex.method.den.toString();
		    checkmethod();
		    status.innerHTML = "";
		},
		error => alert(error)
	    )
	},
	err => alert(err)
    )
}
function checkmethod(){
    let m = glob.simplex.method;
    let checker = (arr, target) => target.every(v => arr.includes(v));
    let ok = checker(m.ions,m.num) & checker(m.ions,m.den)
    togglemethodwarning(ok);
}

function togglemethodwarning(ok){
    let witem = document.getElementById("methodmismatch");
    let iitem = document.getElementById("ions");
    if (ok){
	witem.classList.add("hidden");
	iitem.classList.remove("red");
    } else {
	witem.classList.remove("hidden");
	iitem.classList.add("red");
    }
}

function result2simplex(result){
    glob.simplex = result.data.simplex;
    glob.names = result.data.names;
    glob.class = result.data.class;
    glob.multi = result.data.multi[0];
}

function savejson(){
    let a = document.getElementById('save');
    a.setAttribute("href","data:text/plain," + JSON.stringify(glob));
    a.setAttribute("download","simplex.json");
    a.click();
}

function openjson(){
    let f = document.getElementById("open").files[0];
    let reader = new FileReader();
    reader.onload = function(){
        glob = JSON.parse(reader.result);
	showPresets();
    };
    reader.readAsText(f);
}

// 3. Drift

async function drift(){
    if (glob.start){
	alert("You must run setup first");
	setup();
    } else {
	selectButton(1);
	loadPage("drift.html").then(
	    () => {
		if (glob.class.includes('drift')){ // already drift corrected
		    loadSamples( () => initDrift() )
		} else { // not yet drift corrected
		    loader();
		    shinylight.call("getdrift", {x:glob}, null, extra()).then(
			result => result2simplex(result),
			error => alert(error)
		    ).then(
			() => loadSamples( () => initDrift() ),
			error => alert(error)
		    ).then(
			() => shower(),
			error => alert(error)
		    )
		}
	    },
	    error => alert(error)
	).then(
	    () => loadHelp("help/drift.html"),
	    error => alert(error)
	)
    }
}

function checkratios(){
    glob.ratios = document.getElementById("ratiocheckbox").checked;
}
function checklog(){
    glob.log = document.getElementById("logcheckbox").checked;
}
function checkcov(){
    glob.cov = document.getElementById("covcheckbox").checked;
}
function checkexterr(){
    glob.exterr = document.getElementById("exterrcheckbox").checked;
}
function checkxy(){
    glob.xy = document.getElementById("xycheckbox").checked;
}

function loader(){
    show('.show4loading');
    hide('.hide4loading');
}

function shower(){
    show('.hide4loading');
    hide('.show4loading');
}

function loadSamples(callback){
    let select = document.getElementById("aliquots");
    let samples = glob.simplex.samples;
    let keys = Object.keys(samples);
    for (let i = 0; i<keys.length; i++) {
	let el = document.createElement("option");
	el.textContent = keys[i];
	el.value = i;
	select.appendChild(el);
    }
    callback();
}

function loadTable(dat,header,id,nr){
    let nc = header.length;
    let e = document.getElementById(id);
    e.deg = createDataEntryGrid(id,header,nr);
    e.deg.putCells(0,nr+1,0,nc+1,dat);
}

function initDrift(){
    document.getElementById("aliquots").value = glob.i;
    driftAliquot();
}

function deepcopy(object){
    return(JSON.parse(JSON.stringify(object)))
}

function driftAliquot(){
    glob.i = parseFloat(document.getElementById("aliquots").value);
    let keys = Object.keys(glob.simplex.samples);
    let header = glob.simplex.method.ions;
    let dat = glob.simplex.samples[keys[glob.i]];
    loadTable(dat.time,header,'time-table',dat.time.length);
    loadTable(dat.signal,header,'signal-table',dat.signal.length);
    document.getElementById('outliers').value =
	(dat.outliers===undefined) ? '' : dat.outliers;
}

function backnforth(di,callback){
    let keys = Object.keys(glob.simplex.samples);
    let ns = keys.length;
    glob.i = ((glob.i + di % ns) + ns) % ns; // modulo operation
    document.getElementById("aliquots").value = glob.i;
    callback();
}

function driftPlot(){
    togglehelp(false);
    let keys = Object.keys(glob.simplex.samples);
    let ostring = document.getElementById('outliers').value;
    if (ostring===''){
	delete glob.simplex.samples[keys[glob.i]].outliers;
    } else {
	glob.simplex.samples[keys[glob.i]].outliers = ostring.split(',',10).map(Number);
    }
    shinylight.call('driftPlot', {x:glob},
		    'drift-plot', {'imgType': 'svg'}).then(
			result => {
			    result2simplex(result);
			    shinylight.setElementPlot('drift-plot', result.plot);
			},
			error => alert(error)
		    );
}

function getOutliers(i){
    let e = document.getElementById('drift-plot');
    let omit = e.deg.getColumn(0);
    return(omit);
}

// 3. Logratios

async function logratios(){
    if (glob.start){
	alert("You must run setup first");
	setup();
    } else {
	selectButton(2);
	loadPage("logratios.html").then(
	    () => {
		if (glob.class.includes('logratios')){ // already has logratios
		    loadSamples( () => initLogratios() );
		    document.getElementById("ratiocheckbox").checked = glob.ratios;
		    document.getElementById("logcheckbox").checked = glob.log;
		    if (glob.simplex.method.instrument=='Cameca'){
			show('.show4cameca');
			document.getElementById("xycheckbox").checked = glob.xy;
		    } else {
			hide('.show4cameca');
			glob.xy = false;
		    }
		} else { // does not yet have logratios
		    loader();
		    shinylight.call("getlogratios", {x:glob}, null, extra()).then(
			result => result2simplex(result),
			error => alert(error)
		    ).then(
			() => {
			    loadSamples( () => initLogratios() );
			    document.getElementById("ratiocheckbox").checked = glob.ratios;
			    document.getElementById("logcheckbox").checked = glob.log;
			    document.getElementById("xycheckbox").checked = glob.xy;
			},
			error => alert(error)
		    ).then(
			() => shower(),
			error => alert(error)
		    )
		}
	    },
	    error => alert(error)
	).then(
	    () => loadHelp("help/logratios.html"),
	    error => alert(error)
	)
    }
}

function initLogratios(){
    document.getElementById("aliquots").value = glob.i;
    logratioAliquot();
}

function extra(){
    var extra = {
        'info': function(text) {
            shinylight.setElementText('status', text);
        },
        'progress': function(numerator, denominator) {
            const pc = Math.ceil(numerator * 100 / denominator);
	    let nbars = 20;
	    let ndone = Math.ceil(numerator*nbars/denominator);
	    let done = '|'.repeat(ndone);
	    let todo = '.'.repeat(nbars-ndone);
            shinylight.setElementText('progress', done + todo + ' ' + pc + '%');
        },
	'imgType': 'svg'
    }
    return(extra)
}

function logratioAliquot(){
    glob.i = parseFloat(document.getElementById("aliquots").value);
    let key = Object.keys(glob.simplex.samples)[glob.i];
    let header = glob.names.samples[key].lr.b0g;
    let dat =  glob.simplex.samples[key];
    let b0g = dat.lr.b0g;
    let ns = header.length/2;
    loadTable([b0g.slice(0,ns)],header.slice(0,ns),'b0',1);
    loadTable([b0g.slice(ns,2*ns)],header.slice(ns,2*ns),'g',1);
    document.getElementById('outliers').value =
	(dat.outliers===undefined) ? '' : dat.outliers;
}

function logratioPlot(){
    togglehelp(false);
    show('.plot');
    hide('.table');
    let ostring = document.getElementById('outliers').value;
    let key = Object.keys(glob.simplex.samples)[glob.i];
    if (ostring===''){
	delete glob.simplex.samples[key].outliers;
    } else {
	glob.simplex.samples[key].outliers = ostring.split(',',10).map(Number);
    }
    shinylight.call('logratioPlot',
		    {x:glob, ratios:glob.ratios},
		    'logratio-plot',
		    {'imgType': 'svg'})
	.then(
	    result => shinylight.setElementPlot('logratio-plot', result.plot),
	    error => alert(error)
	);
}

function logratioTable(){
    togglehelp(false);
    show('.table');
    hide('.plot');
    shinylight.call("logratioTable", {x:glob}, null).then(
	result => {
	    let nr = result.data.length;
	    let header = Object.keys(result.data[0]);
	    let tab = createDataEntryGrid('logratio-table', header, nr);
	    shinylight.setGridResult(tab, result);
	},
	error => alert(error)
    );
}

// 4. Calibration

function calibration(){
    if (glob.start){
	alert("You must run setup first");
	setup();
    } else {
	selectButton(3);
	loadPage("calibration.html").then(
	    () => {
		if (glob.class.includes('calibration')){
		    showCalibration();
		} else {
		    createCalibration(showCalibration);
		}
		document.getElementById("shownum").checked = glob.shownum;
	    },
	    error => alert(error)
	).then(
	    () => loadHelp("help/calibration.html"),
	    error => alert(error)
	).then(
	    () => {
		if (glob.calibration.caltype==='regression'){
		    show(".show4geochron");
		    hide(".hide4geochron");
		} else {
		    show(".show4stable");
		    hide(".hide4stable");
		}
	    },
	    error => alert(error)
	)
    }
}
function showCalibration(){
    prepareCalibration();
    // I
    togglecaltype();
    // II
    setstandcomp();
    // III
    setstandsel();
    // set and act upon glob.calibration
}
function prepareCalibration(){
    let cal = glob.calibration;
    if (cal.caltype==null){ // geochron
	cal.caltype = glob.multi ? 'regression' : 'average';
    }
    cal.standtype =
	glob.simplex.calibration.stand.measured[0] ?
	"measured" : "commonradio";
    document.getElementById('caltype').value = cal.caltype;
    document.getElementById('standcomp').value = cal.standcomp;
    document.getElementById('standtype').value = cal.standtype;
    if (cal.preset!==null){
	document.getElementById('presets').value = cal.preset;
    }

}
function createCalibration(callback){
    shinylight.call('createcalibration', {x:glob}, null).then(
	result => result2simplex(result),
	error => alert(error)
    ).then(
	() => callback(),
	error => alert(error)
    );
}
function shownum(){
    glob.shownum = document.getElementById("shownum").checked;
}

// I.

function togglecaltype(){
    let ct = glob.calibration.caltype;
    ct = document.getElementById('caltype').value;
    if (ct === 'average'){
	show('.show4stable');
	hide('.hide4stable');
	if (glob.multi) show('#caltype-warning');
	else hide('#caltype-warning');
    } else { // regression
	show('.show4geochron');
	hide('.hide4geochron');
	setpairing();
	if (glob.multi) hide('#caltype-warning');
	else show('#caltype-warning');
    }
}
function setpairing(){
    let haspairing = glob.simplex.calibration.hasOwnProperty('pairing');
    if (haspairing){
	showpairing();
    } else {
	createPairing(showpairing);
    }
}
function showpairing(){
    let cal = glob.simplex.calibration;
    let nr = cal.pairing.length;
    let header = Object.keys(cal.pairing[0]);
    let val = [new Array(nr)];
    for (let i=0; i<nr; i++){
	val[i] = Object.values(cal.pairing[i]);
    }
    loadTable(val,header,'pairing',nr);
}
function createPairing(callback){
    shinylight.call('createpairing', {x:glob}, null).then(
	result => result2simplex(result),
	error => alert(error)
    ).then(
	() => callback(),
	error => alert(error)
    );
}
function getpairing(){
    let e = document.getElementById('pairing');
    let pairing = glob.simplex.calibration.pairing;
    let nr = pairing.length;
    let header = Object.keys(pairing[0]);
    let dat = e.deg.getCells();
    for (let i=0; i<nr; i++){
	for (let j=0; j<header.length; j++){
	    pairing[i][header[j]] = dat[i][j];
	}
    }
}

// II.
function setstandcomp(){
    let stand = glob.simplex.calibration.stand;
    let header = glob.names.calibration.stand.val;
    let nr = header.length;
    let val = [stand.val];
    let cov = stand.cov;
    togglestandcomp();
    loadTable(val,header,'standlr',1);
    loadTable(cov,header,'standcov',nr);
    togglemismatchwarning();
}
function togglemismatchwarning(){
    let header = glob.names.calibration.stand.val;
    let m = glob.simplex.method;
    let warn = true;
    for (let i=0; i<m.num.length; i++){
	ratio = m.num[i] + '/' + m.den[i];
	if (header.includes(ratio)) {
	    warn = false;
	    break;
	}
    }
    if (warn) show("#mismatch-warning")
    else hide("#mismatch-warning")
}
function togglestandcomp(){
    let sc = glob.calibration.standcomp;
    sc = document.getElementById('standcomp').value;
    if (sc === 'preset2stand'){
	showcalpreset();
	show('.show4preset');
    } else {
	hide('.show4preset')
    }
    if (sc === 't2stand'){
	showtst();
	show('.show4t2stand');
    } else {
	hide('.show4t2stand');
    }
    if (sc === 'd2stand'){
	showcaldel();
	show('.show4d2stand');
    } else {
	hide('.show4d2stand');
    }
}
function showcalpreset(){
    let cal = glob.calibration;
    if (glob.simplex.hasOwnProperty('calibration')){
	let stand = glob.simplex.calibration.stand;
	if (stand.hasOwnProperty('preset')){
	    cal.preset = stand.preset[0];
	    document.getElementById('presets').value = cal.preset;
	}
    }
}
function showtst(){
    let cal = glob.calibration;
    if (glob.simplex.hasOwnProperty('calibration')){
	let stand = glob.simplex.calibration.stand;
	if (stand.hasOwnProperty('tst')){
	    cal.tst[0] = stand.tst[0];
	    cal.tst[1] = stand.tst[1];
	}
    }
    document.getElementById('t').value = cal.tst[0];
    document.getElementById('st').value = cal.tst[1];
}
function showcaldel(){
    let create = true;
    let cal = glob.calibration;
    if (glob.simplex.hasOwnProperty('calibration')){
	let stand = glob.simplex.calibration.stand;
	if (stand.hasOwnProperty('del')){
	    create = false;
	    cal.del.ratios = glob.names.calibration.stand.del.val;
	    cal.del.delval = stand.del.val;
	    cal.del.delcov = stand.del.cov;
	    cal.del.refval = stand.ref.val;
	    cal.del.refcov = stand.ref.cov;
	}
	if (stand.hasOwnProperty('ref') && stand.ref.hasOwnProperty('preset')){
	    document.getElementById('deltaref').value = stand.ref.preset[0];
	}
    }
    let elr = document.getElementById('standlr');
    if (create){
	let header = elr.deg.getColumnHeaders();
	cal.del.ratios = header;
	let edel = document.getElementById('deltab');
	let ecov = document.getElementById('delcovtab');
	let eref = document.getElementById('delreftab');
	edel.deg = createDataEntryGrid('deltab',header,1);
	ecov.deg = createDataEntryGrid('delcovtab',header,header.length);
	eref.deg = createDataEntryGrid('delreftab',header,1);
    } else {
	let header = cal.del.ratios;
	loadTable([cal.del.delval],header,'deltab',1);
	loadTable(cal.del.delcov,header,'delcovtab',header.length);
	loadTable([cal.del.refval],header,'delreftab',1);
    }
}
function togglestandtype(){
    glob.calibration.standtype = document.getElementById('standtype').value;
    createCalibration(showCalibration);
}
function preset2standard(){
    glob.calibration.preset = document.getElementById('presets').value;
    shinylight.call('preset2standard', {x:glob}, null).then(
	result => {
	    result2simplex(result);
	    setstandcomp();
	},
	error => alert(error)
    )
}
function t2stand(){
    let cal = glob.calibration;
    cal.tst[0] = parseFloat(document.getElementById('t').value);
    cal.tst[1] = parseFloat(document.getElementById('st').value);
    shinylight.call('t2stand', {x:glob}, null).then(
	result => {
	    result2simplex(result);
	    setstandcomp();
	},
	error => alert(error)
    )
}
function d2stand(){
    let cal = glob.calibration;
    let edel = document.getElementById('deltab');
    let ecov = document.getElementById('delcovtab');
    let eref = document.getElementById('delreftab');
    cal.del.delval = edel.deg.getCells();
    cal.del.delcov = ecov.deg.getColumns();
    cal.del.refval = eref.deg.getCells();
    shinylight.call('d2stand', {x:glob}, null).then(
	result => {
	    result2simplex(result);
	    setstandcomp();
	},
	error => alert(error)
    )
}
function deltaref(){
    let ref = document.getElementById('deltaref').value;
    shinylight.call('deltaref', {ref:ref}, null).then(
	result => {
	    result2simplex(result);
	    setstandcomp();
	},
	error => alert(error)
    )
}

// III.
function setstandsel(){
    let cal = glob.simplex.calibration;
    let hasprefix = cal.hasOwnProperty('prefix');
    if (!hasprefix) cal.prefix = '';
    if (glob.standards.length<1){
	prefix2standards();
    }
    document.getElementById('prefix').value = cal.prefix;
    markStandards();
}
function prefix2standards(){
    let keys = Object.keys(glob.simplex.samples);
    let prefix = glob.simplex.calibration.prefix;
    glob.standards = [];
    for (let i=0; i<keys.length; i++){
	if (keys[i].indexOf(prefix) !== -1){
	    glob.standards.push(keys[i]);
	}
    }
}
function updateStandardPrefix(){
    glob.simplex.calibration.prefix = document.getElementById('prefix').value;
    prefix2standards();
    markStandards();
}
function markStandards(){
    let keys = Object.keys(glob.simplex.samples);
    let nk = keys.length;
    let dat = new Array(nk);
    for (let i=0; i<nk; i++){
	if (glob.standards.includes(keys[i])){
	    dat[i] = [keys[i],'yes'];
	} else {
	    dat[i] = [keys[i],'no'];
	}
    }
    loadTable(dat,['aliquots','selected?'],'aliquots',nk);
}

// IV.
function calibrator(){
    togglehelp(false);
    registerStandards();
    shinylight.call('calibrator', {x:glob},
		    'calibration-plot', {'imgType': 'svg'}).then(
	result => {
	    result2simplex(result),
	    shinylight.setElementPlot('calibration-plot', result.plot)
	},
	error => alert(error)
    )
}
function registerStandards(){
    let stand = glob.simplex.calibration.stand;
    let e = document.getElementById('standlr');
    stand.val = e.deg.getCells()[0].map(Number);
    e = document.getElementById('standcov');
    let cov = e.deg.getCells();
    for (let i=0; i<stand.val.length; i++){
	stand.cov[i] = cov[i].map(Number);
    }
    if (glob.calibration.caltype=='regression'){
	e = document.getElementById('pairing');
	let pairing = glob.simplex.calibration.pairing;
	let content = e.deg.getCells();
	let header = e.deg.getColumnHeaders();
	for (let i=0; i<header.length; i++){
	    pairing[0][header[i]] = content[0][i];
	}
    }
    e = document.getElementById('aliquots');
    let dat = e.deg.getColumns();
    glob.standards = [];
    for (let i=0; i<dat.aliquots.length; i++){
	if (dat['selected?'][i]==='yes') glob.standards.push(dat.aliquots[i]);
    }
}

// 5. samples

function samples(){
    if (glob.start){
	alert("You must run setup first");
	setup();
    } else {
	selectButton(4);
	loadPage("samples.html").then(
	    () => {
		setsampsel();
		document.getElementById("shownum").checked = glob.shownum;
		document.getElementById("logcheckbox").checked = glob.log;
		document.getElementById("covcheckbox").checked = glob.cov;
		document.getElementById("exterrcheckbox").checked = glob.exterr;
	    }, error => alert(error)
	).then(
	    () => loadHelp("help/samples.html"),
	    error => alert(error)
	)
    }
}
function setsampsel(){
    document.getElementById('prefix').value = glob.sampleprefix.join(',');
    if (glob.samples.length<1){
	prefix2samples();	
    }
    markSamples();
}
function prefix2samples(){
    let keys = Object.keys(glob.simplex.samples);
    glob.samples = [];
    for (let i=0; i<keys.length; i++){
	for (let j=0; j<glob.sampleprefix.length; j++){
	    if (keys[i].indexOf(glob.sampleprefix[j]) !== -1){
		glob.samples.push(keys[i]);
	    }
	}
    }
}
function updateSamplePrefix(){
    let prefixes = document.getElementById('prefix').value;
    glob.sampleprefix = prefixes.split(',');
    prefix2samples();
    markSamples();
}
function markSamples(){
    let keys = Object.keys(glob.simplex.samples);
    let nk = keys.length;
    let dat = new Array(nk);
    for (let i=0; i<nk; i++){
	if (glob.samples.includes(keys[i])){
	    dat[i] = [keys[i],'yes'];
	} else {
	    dat[i] = [keys[i],'no'];
	}
    }
    loadTable(dat,['aliquots','selected?'],'aliquots',nk);
}
function calibrate(plot){
    togglehelp(false);
    registerSamples();
    if (plot) calibrate_plot()
    else calibrate_table()
}
function registerSamples(){
    let e = document.getElementById('aliquots');
    let dat = e.deg.getColumns();
    glob.samples = [];
    for (let i=0; i<dat.aliquots.length; i++){
	if (dat['selected?'][i]==='yes') glob.samples.push(dat.aliquots[i]);
    }
}
function calibrate_plot(){
    show('.plot');
    hide('.table');
    shinylight.call("calibrateSamples",{x:glob},'sample-calibration-plot',{'imgType': 'svg'})
	.then(
	    result => shinylight.setElementPlot('sample-calibration-plot', result.plot),
	    error => alert(error)
	)
}
function calibrate_table(){
    show('.table');
    hide('.plot');
    shinylight.call("calibratedTable", {x:glob}, null).then(
	result => {
	    let nr = result.data.length;
	    let header = Object.keys(result.data[0]);
	    let tab = createDataEntryGrid('sample-calibration-table', header, nr);
	    shinylight.setGridResult(tab, result);
	},
	error => alert(error)
    );
}

// 6. finish
function finish(){
    if (glob.start){
	alert("You must run setup first");
	setup();
    } else {
	selectButton(5);
	loadPage("finish.html").then(
	    () => finishinit(),
	    error => alert(error)
	).then(
	    () => loadHelp("help/finish.html"),
	    error => alert(error)
	).then(
	    () => {
		if (glob.calibration.caltype==='regression'){
		    show(".show4geochron");
		    hide(".hide4geochron");
		} else {
		    show(".show4stable");
		    hide(".hide4stable");
		}
	    },
	    error => alert(error)
	)
    }
}
function finishinit(){
    let cal = glob.calibration;
    if (cal.caltype=='average'){
	show('.show4stable');
	hide('.hide4stable');
	document.getElementById('deltatype').value = glob.delta.type;
	showdeltasettings();
    } else {
	show('.show4geochron');
	hide('.hide4geochron');
	document.getElementById('IsoplotRtype').value = glob.IsoplotRtype;
    }
}

function showdeltasettings(){
    let del = glob.delta;
    if (del.preset==null & del.ratios==null & del.val==null){
	let stand = glob.simplex.calibration.stand;
	if (stand.ref.hasOwnProperty('preset')){
	    del.preset = stand.ref.preset[0];
	    document.getElementById('deltaref').value = del.preset;
	}
	del.ratios = glob.names.calibration.stand.ref.val.slice();
	del.val = stand.ref.val.slice();
    }
    if (['VSMOW','troilite'].includes(del.preset)){
	document.getElementById('deltaref').value = del.preset;
    }
    loadTable([del.val],del.ratios,'delreftab',1);
}
function preset2deltaref(){
    glob.delta.preset = document.getElementById('deltaref').value;
    shinylight.call('preset2deltaref', {x:glob}, null).then(
	result => loadTable([result.data.val],result.data.ratios,'delreftab',1),
	error => alert(error)
    )
}

function convert(fn){
    togglehelp(false);
    if (glob.calibration.caltype=='average'){
	let e = document.getElementById('delreftab');
	glob.delta.val = e.deg.getCells()[0].map(Number);
    }
    shinylight.call(fn, {x:glob}, null).then(
	result => {
	    let nr = result.data.length;
	    let header = Object.keys(result.data[0]);
	    let tab = createDataEntryGrid('final-table', header, nr);
	    shinylight.setGridResult(tab, result);
	},
	error => alert(error)
    );
}

function toggledeltatype(){
    glob.delta.type = document.getElementById('deltatype').value;
}

function toggleIsoplotRtype(){
    glob.IsoplotRtype = document.getElementById('IsoplotRtype').value;
}

function export2IsoplotR(){
    glob.IsoplotRformat = document.getElementById("IsoplotRtype").value;
    var json = null;
    fetch('js/IsoplotR.json')
	.then(response => {
	    if (!response.ok) {
		throw new Error("HTTP error " + response.status);
	    } else {
		return response.json();
	    }
	}).then(
	    result => json = result,
	    err => alert(err)
	).then(
	    async () => {
		let result = await shinylight.call('export2IsoplotR',
						   { x:glob }, null);
		let gc = null;
		let pd = null;
		let format = null;
		switch (glob.IsoplotRtype){
		case 'U-Pb':
		    gc = 'U-Pb';
		    pd = 'concordia';
		    format = 5;
		    break;
		case 'U-Th-Pb':
		    gc = 'U-Pb';
		    pd = 'concordia';
		    format = 8;
		    break;		    
		case 'Th-Pb':
		    gc = 'Th-Pb';
		    pd = 'isochron';
		    format = 2;
		    break;		    
		}
		json.settings.geochronometer = gc;
		json.settings.plotdevice = pd;
		json.settings[gc].format = format;
		json.data[gc] = result;
		let nr = result.data[Object.keys(result.data)[0]].length;
		json.data[gc].data['(C)'] = Array.from(' '.repeat(nr));
		json.data[gc].data['(omit)'] = Array.from(' '.repeat(nr));
		json.data[gc].ierr = 1;
	    },
	    err => alert(err)
	)
	.then(
	    () => {
		let fname = prompt("Please enter a file name", "IsoplotR.json");
		if (fname != null){
		    document.getElementById('fname').setAttribute(
			"href","data:text/plain," + JSON.stringify(json)
		    );
		    document.getElementById('fname').setAttribute("download",fname);
		    document.getElementById('fname').click();
		}
	    },
	    err => alert(err)
	);
}
