// 1. Startup

var simplex = {
    'selected': 0,
    'buttonIDs': ['setup','drift','logratios','calibration','samples','finish'],
    'method': null,
    'samples': null
}

var elements = {
    'driftime': null,
    'driftsig': null
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
        simplex.buttonIDs[simplex.selected]
    );
}

function selectButton(i){
    selectedButton().classList.remove('on')
    simplex.selected = i;
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
        () => loadPresets()
    ).catch(
        error => alert(error)
    );
}

async function loadPresets(){
    const m = document.getElementById("methods").value;
    const result = await shinylight.call('presets', { method: m }, null);
    simplex.samples = result.data.samples;
    simplex.method = result.data.method;
    showPresets();
    fileFormats();
}

function showPresets(){
    let assign = (id) => document.getElementById(id).value = simplex.method[id];
    assign('description');
    assign('instrument');
    assign('ions');
    assign('num');
    assign('den');
    let stabset = document.querySelector('.hide4stable');
    if (stable()) {
	stabset.classList.add('hidden')
    } else {
	stabset.classList.remove('hidden')
	assign('oxide');
    }
    document.getElementById('nominalblank').checked =
	simplex.method.nominalblank[0];
}

function stable(){
    let m = simplex.method.method;
    return(["IGG-O","IGG-S","GA-O"].includes(m[0]))
}

function fileFormats(){
    let accept = ['.asc','.op','.pd'];
    if (simplex.method.instrument=='Cameca'){
	accept = '.asc';
    } else if (simplex.method.instrument=='SHRIMP'){
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
		result => simplex.samples = result.data.samples,
		error => alert(error)
	    )
	},
	err => alert(err)
    )
}

// 3. Drift

function drift(){
    selectButton(1);
    loadPage("drift.html").then(
	() => loadSamples(),
	error => alert(error)
    ).then(
	shinylight.call("driftCorr", {x:simplex}, null).then(
	    result => simplex.samples = result.data.samples,
	    error => alert(error)
	)
    );
}

function loadSamples(){
    let select = document.getElementById("aliquots");
    let samples = simplex.samples;
    let keys = Object.keys(samples);
    for(let i = 0; i<keys.length; i++) {
	let el = document.createElement("option");
	el.textContent = keys[i];
	el.value = i;
	select.appendChild(el);
    }
    let nr = samples[keys[0]].signal.length;
    let header = simplex.method.ions;
    elements.driftime = createDataEntryGrid('drift-time-table',header, nr);
    elements.driftsig = createDataEntryGrid('drift-signal-table',header,nr);
    driftTable(0);
}

function driftTable(i){
    let samples = simplex.samples;
    let keys = Object.keys(samples);
    let nr = elements.driftime.rowCount();
    let nc = elements.driftsig.columnCount();
    let tim = stringtab(samples[keys[i]].time);
    let sig = stringtab(samples[keys[i]].signal);
    elements.driftime.putCells(0,nr+1,0,nc+1,tim);
    elements.driftsig.putCells(0,nr+1,0,nc+1,sig);
}

function stringtab(dat){
    let nr = dat.length;
    let nc = dat[0].length;
    let arr = new Array(nr);
    for (let i=0; i<nr; i++){
	arr[i] = new Array(nc);
	for (let j=0; j<nc; j++){
	    arr[i][j] = parseFloat(dat[i][j]);
	}
    }
    return(arr)
}

function showAliquot(){
    let aliquot = document.getElementById("aliquots").value;
    driftTable(aliquot);
}

function driftPlot(){
    let i = document.getElementById("aliquots").value;
    shinylight.call('driftPlot', {x:simplex, i:i}, 'drift-plot').then(
	result => shinylight.setElementPlot('drift-plot', result.plot),
	error => alert(error)
    );
}

// 3. Logratios

function logratios(){
    selectButton(2);
    loadPage("logratios.html");
}

function calibration(){
    selectButton(3);
    loadPage("calibration.html");
}

function samples(){
    selectButton(4);
    loadPage("samples.html");
}

function finish(){
    selectButton(5);
    loadPage("finish.html");
}
