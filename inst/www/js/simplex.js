// 1. Startup

var glob = {
    'simplex': {
	'method': null,
	'samples': null,
	'standard': null,
	'calibration': null,
	'calibrated': null
    },
    'names': null,
    'selected': 0,
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
        () => loadPresets(),
        error => alert(error)
    );
}

async function loadPresets(){
    const m = document.getElementById("methods").value;
    const result = await shinylight.call('presets', { method: m }, null);
    result2simplex(result);
    showPresets();
    fileFormats();
}

function showPresets(){
    let assign = (id) => document.getElementById(id).value =
	glob.simplex.method[id];
    assign('description');
    assign('instrument');
    assign('ions');
    assign('num');
    assign('den');
    showOrHide('.hide4stable',stable(),assign,'oxide');
    showOrHide('.hide4multi',glob.simplex.method.multicollector[0]);
    labelButtons();
    document.getElementById('multicollector').checked =
	glob.simplex.method.multicollector[0];
    document.getElementById('nominalblank').checked =
	glob.simplex.method.nominalblank[0];
}

function showOrHide(cls,condition,callback,arg){
    let set = document.querySelector(cls);
    if (condition) {
	set.classList.add('hidden');
    } else {
	set.classList.remove('hidden');
	if (callback !== undefined && arg !== undefined) callback(arg)
    }
}

function stable(){
    let m = glob.simplex.method.method;
    return(["IGG-O","IGG-S","GA-O"].includes(m[0]))
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
    )
}

function result2simplex(result){
    glob.simplex = result.data.simplex;
    glob.names = result.data.names;
}

// 3. Drift

function drift(){
    selectButton(1);
    loadPage("drift.html").then(
	() => loader(),
	error => alert(error)
    ).then(
	shinylight.call("getdrift", {x:glob}, null).then(
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
    shinylight.call('driftPlot', {x:glob, i:i}, 'drift-plot').then(
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
	shinylight.call("getlogratios", {x:glob}, null).then(
	    result => result2simplex(result),
	    error => alert(error)
	).then(
	    () => loadSamples(() => logratioAliquot()),
	    error => alert(error)
	).then(
	    () => shower(),
	    error => alert(error)
	)
    )
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
    shinylight.call('logratioPlot', {x:glob, i:i}, 'logratio-plot').then(
	result => shinylight.setElementPlot('logratio-plot', result.plot),
	error => alert(error)
    );
}

// 4. Calibration

function calibration(){
    selectButton(3);
    loadPage("calibration.html").then(
	() => loader(),
	error => alert(error)
    ).then(
	() => {
	    showStandard();
	    showPrefix();
	},
	error => alert(error)
    ).then(
	() => showPrefixed(),
	error => alert(error)
    ).then(
	() => shower(),
	error => alert(error)
    );
}

function showStandard(){
    document.getElementById('standards').value = glob.simplex.standard.name[0];
}

function showPrefix(){
    document.getElementById('prefix').value = glob.simplex.standard.prefix;
}

function showPrefixed(){
    let keys = Object.keys(glob.simplex.samples);
    let dat = new Array(keys.length);
    for (let i=0; i<keys.length; i++){
	dat[i] = [keys[i],'yes'];
    }
    loadTable(dat,['aliquots','selected?'],'aliquots',keys.length);
}

function chooseStandard(){
    let i = document.getElementById("standards").value;
    let key = Object.keys(glob.simplex.samples)[i];
    loadTable(null,'','aliquots',1);
    loadTable(dat,header,id,nr)
}

function samples(){
    selectButton(4);
    loadPage("samples.html");
}

function finish(){
    selectButton(5);
    loadPage("finish.html");
}
