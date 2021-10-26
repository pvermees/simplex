// 1. Startup

var simplex = {
    'selected': 0,
    'buttonIDs': ['setup','drift','logratios','calibration','samples','finish'],
    'method': null,
    'names': null,
    'samples': null,
    'logratios': false
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
    result2simplex(result);
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
    showOrHide('.hide4stable',stable(),assign,'oxide');
    showOrHide('.hide4multi',simplex.method.multicollector[0]);
    labelButtons();
    document.getElementById('multicollector').checked =
	simplex.method.multicollector[0];
    document.getElementById('nominalblank').checked =
	simplex.method.nominalblank[0];
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
    let m = simplex.method.method;
    return(["IGG-O","IGG-S","GA-O"].includes(m[0]))
}

function labelButtons(){
    let labelNums = simplex.method.multicollector[0] ?
	[1,0,2,3,4,5] : [1,2,3,4,5,6];
    let id = null;
    for (let i=0; i<labelNums.length; i++){
	id = simplex.buttonIDs[i];
	document.getElementById(id).innerText =
	    labelNums[i] + '. ' + simplex.buttonIDs[i];
    }
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
		result => result2simplex(result),
		error => alert(error)
	    )
	},
	err => alert(err)
    )
}

function result2simplex(result){
    simplex.samples = result.data.samples;
    simplex.names = result.data.names;
    simplex.method = result.data.method;
}

// 3. Drift

function drift(){
    selectButton(1);
    loadPage("drift.html").then(
	() => loader(),
	error => alert(error)
    ).then(
	shinylight.call("getdrift", {x:simplex}, null).then(
	    result => result2simplex(result),
	    error => alert(error)
	).then(
	    () => loadDriftSamples(),
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

function loadDriftSamples(){
    let select = document.getElementById("aliquots");
    let samples = simplex.samples;
    let keys = Object.keys(samples);
    for (let i = 0; i<keys.length; i++) {
	let el = document.createElement("option");
	el.textContent = keys[i];
	el.value = i;
	select.appendChild(el);
    }
    driftAliquot();
}

function loadTable(dat,header,id){
    let nr = dat.length;
    let nc = header.length;
    let tab = createDataEntryGrid(id,header,nr);
    tab.putCells(0,nr+1,0,nc+1,dat);
}

function driftAliquot(){
    let i = document.getElementById("aliquots").value;
    let keys = Object.keys(simplex.samples);
    let header = simplex.method.ions;
    loadTable(simplex.samples[keys[i]].time,header,'time-table');
    loadTable(simplex.samples[keys[i]].signal,header,'signal-table');
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
    loadPage("logratios.html").then(
	() => loader(),
	error => alert(error)
    ).then(
	shinylight.call("getlogratios", {x:simplex}, null).then(
	    result => result2simplex(result),
	    error => alert(error)
	).then(
	    () => loadLogratioSamples(),
	    error => alert(error)
	).then(
	    () => shower(),
	    error => alert(error)
	).then(
	    () => document.getElementById("logratios").checked =
		simplex.logratios
	)
    )
}

function checkLogratios(){
    simplex.logratios = document.getElementById("logratiocheckbox").checked;
}

function loadLogratioSamples(){
    let select = document.getElementById("aliquots");
    let samples = simplex.samples;
    let keys = Object.keys(samples);
    for(let i = 0; i<keys.length; i++) {
	let el = document.createElement("option");
	el.textContent = keys[i];
	el.value = i;
	select.appendChild(el);
    }
    logratioAliquot();
}

function logratioAliquot(){
    let i = document.getElementById("aliquots").value;
    let keys = Object.keys(simplex.samples);
    let header = simplex.method.num;
    loadTable(simplex.samples[keys[i]].lr.t,header,'time-table');
    loadTable(simplex.samples[keys[i]].lr.obs,header,'signal-table');
}

function logratioPlot(){
    let i = document.getElementById("aliquots").value;
    shinylight.call('logratioPlot',
		    {x:simplex, i:i, logratios:simplex.logratios},
		    'logratio-plot').then(
	result => shinylight.setElementPlot('logratio-plot', result.plot),
	error => alert(error)
    );
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
