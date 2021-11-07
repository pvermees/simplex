// 1. Startup

var glob = {
    'simplex': {
	'method': { 'method': ['IGG-UPb'] },
	'samples': null,
	'standard': null,
	'calibration': null,
	'calibrated': null
    },
    'names': null,
    'class': 'simplex',
    'selected': 0,
    'ratios': false,
    'sampleprefix': null,
    'standards': [],
    'datatype': null,
    'IsoplotRformat': 'U-Pb',
    'buttonIDs': ['setup','drift','logratios','calibration','samples','finish']
}

function start() {
    // This will be moved into shinylight.initialize()
    new Promise((resolve, reject) => {
        rrpc.initialize(() => resolve(), error => reject(error));
    }).then(() =>
        setup()
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

async function loadPage(url) {
    let response = await fetch(url);
    let text = await response.text();
    document.getElementById("contents").innerHTML = text;
}

// 2. Setup

function setup(){
    selectButton(0);
    loadPage("setup.html").then(
        () => {
	    initpreset();
	    loadPresets()
	},
        error => alert(error)
    );
}

function initpreset(){
    document.getElementById("methods").value = glob.simplex.method.method[0];
}

async function loadPresets(){
    const m = document.getElementById("methods").value;
    glob.simplex.method.method[0] = m;
    const result = await shinylight.call('presets', { method: m }, null);
    result2simplex(result);
    showPresets();
    fileFormats();
    setDataType();
}

function setDataType(){
    shinylight.call('getdatatype', { x: glob }, null).then(
	result => glob.datatype = result.data[0],
	error => alert(error)
    )
}

function showPresets(){
    let assign = (id) => document.getElementById(id).value =
	glob.simplex.method[id];
    assign('description');
    assign('instrument');
    assign('ions');
    assign('num');
    assign('den');
    hideIt('.hide4stable',stable(),assign,'oxide');
    hideIt('.hide4multi',glob.simplex.method.multicollector[0]);
    labelButtons();
    document.getElementById('multicollector').checked =
	glob.simplex.method.multicollector[0];
    document.getElementById('nominalblank').checked =
	glob.simplex.method.nominalblank[0];
}

function showIt(cls,condition,callback,arg){
    showOrHide(cls,condition,callback,arg)
}
function hideIt(cls,condition,callback,arg){
    showOrHide(cls,!condition,callback,arg)
}
function showOrHide(cls,condition,callback,arg){
    let set = document.querySelector(cls);
    if (set==null) return
    if (condition){
	set.classList.remove('hidden');
	if (callback !== undefined && arg !== undefined) callback(arg)
    } else {
	set.classList.add('hidden');
    }
}

function stable(){
    return(["oxygen","sulphur"].includes(glob.datatype))
}

function geochron(){
    return(["U-Pb","U-Th-Pb"].includes(glob.datatype))
}

function labelButtons(){
    let labelNums = glob.simplex.method.multicollector[0] ?
	[1,0,2,3,4,5] : [1,2,3,4,5,6];
    let id = null;
    for (let i=0; i<labelNums.length; i++){
	id = glob.buttonIDs[i];
	document.getElementById(id).innerText =
	    labelNums[i] + '. ' + glob.buttonIDs[i];
    }
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
    let status = document.getElementById('upload-status');
    let f = document.getElementById('upload').files;
    let fns = {};
    let tcs = {};
    for (let i=0; i<f.length; i++){
	status.innerHTML = " Loading file " + (i+1) + " of " + f.length;
	fns[i] = f[i].name;
	tcs[i] = await readFile(f[i]);
	status.innerHTML = (i==f.length-1) ? "" :
	    " Loaded file " + (i+1) + " of " + f.length;

    }
    return({fns:fns, tcs:tcs})
}

async function upload(){
    const m = document.getElementById("methods").value;
    readFiles().then(
	f => {
	    shinylight.call('upload', {f:f, m:m}, null).then(
		result => result2simplex(result),
		error => alert(error)
	    )
	},
	err => alert(err)
    ).then(
	() => setDataType(),
	err => alert(err)
    )
}

function result2simplex(result){
    glob.simplex = result.data.simplex;
    glob.names = result.data.names;
    glob.class = result.data.class;
}

// 3. Drift

function drift(){
    selectButton(1);
    loadPage("drift.html").then(
	() => loader(),
	error => alert(error)
    ).then(
	shinylight.call("getdrift", {x:glob}, null, extra()).then(
	    result => result2simplex(result),
	    error => alert(error)
	).then(
	    () => loadSamples(() => driftAliquot()),
	    error => alert(error)
	).then(
	    () => shower(),
	    error => alert(error)
	)
    )
}

function checkratios(){
    glob.ratios = document.getElementById("ratiocheckbox").checked;
}

function loader(){
    let show = document.querySelector('.show4loading');
    let hide = document.querySelector('.hide4loading');
    show.classList.remove('hidden');
    hide.classList.add('hidden');
}

function shower(){
    let show = document.querySelector('.show4loading');
    let hide = document.querySelector('.hide4loading');
    show.classList.add('hidden');
    hide.classList.remove('hidden');
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
    let tab = createDataEntryGrid(id,header,nr);
    tab.putCells(0,nr+1,0,nc+1,dat);
}

function driftAliquot(){
    let i = document.getElementById("aliquots").value;
    let keys = Object.keys(glob.simplex.samples);
    let header = glob.simplex.method.ions;
    let dat = glob.simplex.samples[keys[i]];
    loadTable(dat.time,header,'time-table',dat.time.length);
    loadTable(dat.signal,header,'signal-table',dat.signal.length);
}

function driftPlot(){
    let i = document.getElementById("aliquots").value;
    shinylight.call('driftPlot',
		    {x:glob, i:i},
		    'drift-plot',
		    {'imgType': 'svg'}).then(
	result => shinylight.setElementPlot('drift-plot', result.plot),
	error => alert(error)
    );
}

// 3. Logratios

function logratios(){
    selectButton(2);
    loadPage("logratios.html").then(
	() => loader(),
	error => alert(error)
    ).then(
	shinylight.call("getlogratios", {x:glob}, null, extra()).then(
	    result => result2simplex(result),
	    error => alert(error)
	).then(
	    () => {
		loadSamples(() => logratioAliquot());
		document.getElementById("ratiocheckbox").checked = glob.ratios;
	    },
	    error => alert(error)
	).then(
	    () => shower(),
	    error => alert(error)
	)
    )
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
    let i = document.getElementById("aliquots").value;
    let key = Object.keys(glob.simplex.samples)[i];
    let header = glob.names.samples[key].lr.b0g;
    let b0g = glob.simplex.samples[key].lr.b0g;
    let ns = header.length/2;
    loadTable([b0g.slice(0,ns)],header.slice(0,ns),'b0',1);
    loadTable([b0g.slice(ns,2*ns)],header.slice(ns,2*ns),'g',1);
}

function logratioPlot(){
    let i = document.getElementById("aliquots").value;
    shinylight.call('logratioPlot',
		    {x:glob, i:i, ratios:glob.ratios},
		    'logratio-plot',
		    {'imgType': 'svg'})
	.then(
	    result => shinylight.setElementPlot('logratio-plot', result.plot),
	    error => alert(error)
	);
}

// 4. Calibration

function calibration(){
    selectButton(3);
    loadPage("calibration.html").then(
	() => {
	    showOrHideStandards();
	    if (typeof glob.simplex.standard != 'undefined'){
		document.getElementById('standards').value =
		    glob.simplex.standard.name[0];
		document.getElementById('prefix').value =
		    glob.simplex.standard.prefix;
	    } else {
		glob.simplex.standard = {
		    'name': [''],
		    'prefix': ''
		}
	    }
	    markStandardsByPrefix()
	},
	error => alert(error)
    );
}

function showOrHideStandards(){
    let disable = null;
    switch(glob.datatype) {
    case 'U-Pb':
	disable = [5,6];
        break; 
    case 'U-Th-Pb':
	disable = [1,2,5,6];
	break;
    case 'oxygen': 
        disable = [1,2,3,4,6];
        break;
    case 'sulphur':
	disable = [1,2,3,4,5];
        break; 
    default: //optional
	//statements
    }
    for (let i=0; i<disable.length; i++){
	document.getElementById("standards").options[disable[i]].disabled = true;
    }
}

function markStandardsByPrefix(){
    let prefix = document.getElementById('prefix').value;
    glob.simplex.standard.prefix = prefix;
    let keys = Object.keys(glob.simplex.samples);
    let nk = keys.length;
    let dat = new Array(nk);
    glob.standards = new Array(nk);
    for (let i=0; i<nk; i++){
	if (keys[i].indexOf(prefix) !== -1){
	    glob.standards[i] = keys[i];
	    dat[i] = [keys[i],'yes'];
	} else {
	    dat[i] = [keys[i],'no'];
	}
    }
    loadTable(dat,['aliquots','selected?'],'aliquots',keys.length);
}

function chooseStandard(){
    let stand = document.getElementById("standards").value;
    shinylight.call("getstandard", {preset:stand}, null).then(
	result => glob.simplex.standard = result.data,
	error => alert(error)
    )
}

function calibrator(){
    shinylight.call("calibrator", {x:glob},
		    'calibration-plot', {'imgType': 'svg'}).then(
	result => {
	    result2simplex(result),
	    shinylight.setElementPlot('calibration-plot', result.plot)
	},
	error => alert(error)
    )
}

// 5. samples

function samples(){
    selectButton(4);
    loadPage("samples.html").then(
	() => markSamplesByPrefix(),
	error => alert(error)
    );
}

function markSamplesByPrefix(){
    glob.sampleprefix = document.getElementById('prefix').value;
    let keys = Object.keys(glob.simplex.samples);
    let nk = keys.length;
    let dat = new Array(nk);
    for (let i=0; i<nk; i++){
	if (keys[i].indexOf(glob.sampleprefix) !== -1){
	    dat[i] = [keys[i],'yes'];
	} else {
	    dat[i] = [keys[i],'no'];
	}
    }
    loadTable(dat,['aliquots','selected?'],'aliquots',keys.length);
}

function calibrate(){
    shinylight.call("calibrateSamples",
		    {x:glob},
		    'sample-calibration-plot',
		    {'imgType': 'svg'}).then(
	result => shinylight.setElementPlot('sample-calibration-plot', result.plot),
	error => alert(error)
    )
}

// 6. finish

function finish(){
    selectButton(5);
    loadPage("finish.html").then(
	() => {
	    hideIt('.hide4stable',stable());
	    showIt('.show4UThPb',glob.datatype==='U-Th-Pb');
	    document.getElementById('prefix').value = glob.sampleprefix;
	    markSamplesByPrefix();
	}, error => alert(error)
    );
}

function plotresults(){
    hideIt('.hide4plot',true);
    hideIt('.hide4table',false);
    shinylight.call("plotresults",
		    {x:glob},
		    'final-plot',
		    {'imgType': 'svg'}).then(
	result => {
	    shinylight.setElementPlot('final-plot', result.plot)
	},
	error => alert(error)
    );
}

function resultstable(){
    hideIt('.hide4table',true);
    hideIt('.hide4plot',false);
    shinylight.call("resultstable", {x:glob}, null).then(
	result => {
	    let nr = result.data.length;
	    let header = Object.keys(result.data[0]);
	    let tab = createDataEntryGrid('final-table', header, nr);
	    shinylight.setGridResult(tab, result);
	},
	error => alert(error)
    );
}

function export2isoplotr(){
    glob.IsoplotRformat = document.getElementById("format").value;
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
		let result = await shinylight.call('export2isoplotr',
						   { x:glob }, null);
		let gc = null;
		let pd = null;
		let format = null;
		switch (glob.IsoplotRformat){
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
	    },
	    err => alert(err)
	)
	.then(
	    () => {
		let fname = prompt("Please enter a file name", "simplex.json");
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
