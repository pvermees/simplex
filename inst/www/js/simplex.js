var simplex = {
    'selected': 0,
    'buttonIDs': ['setup','drift','logratios','calibration','samples','finish'],
    'method': null,
    'data': null
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
    simplex.method = result.data;
    showPresets();
}

function showPresets(){
    let assign = (id) => document.getElementById(id).value = simplex.method[id];
    assign('description');
    assign('instrument');
    assign('ions');
    assign('num');
    assign('den');
    document.getElementById('nominalblank').checked = simplex.method.nominalblank[0];
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

async function upload(){
    let f = document.getElementById('upload').files;
    let status = document.getElementById('upload-status');
    let txt = null;
    const m = document.getElementById("methods").value;
    for (let i=0; i<f.length; i++){
	status.innerHTML = " Loading file " + (i+1) + " of " + f.length;
	txt = await readFile(f[i]);
	shinylight.call('upload', {f: txt, m: m}, null).then(
	    result => {
		simplex.data = result.data;
		status.innerHTML = (i==f.length-1) ? "" :
		    " Loaded file " + (i+1) + " of " + f.length;
	    },
	    error => alert(error)
	);
    }
}

function drift(){
    selectButton(1);
    loadPage("drift.html");
}

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
