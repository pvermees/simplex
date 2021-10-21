var simplex = {
    'selected': 0,
    'buttonIDs': ['setup','drift','logratios','calibration','samples','finish'],
    'oncol': 'yellow',
    'offcol': 'white',
    'method': null,
    'data': null
}

function start() {
    rrpc.initialize();
    shinylight.initialize();
    setup();
}

function selectButton(i){
    document.getElementById(simplex.buttonIDs[simplex.selected]).style.background =
	simplex.offcol;
    document.getElementById(simplex.buttonIDs[i]).style.background =
	simplex.oncol;
    simplex.selected = i;
}

function loadPage(url) {
    fetch(url)
	.then(response => response.text())
	.then(text => {
	    document.getElementById("contents").innerHTML = text;
	});
}

function setup(){
    selectButton(0);
    loadPage("setup.html"); // TODO: load presets when page has loaded
}

function loadPresets(){
    let m = document.getElementById("methods").value;
    shinylight.call('presets', {method: m}, null).then(
	result => {
	    simplex.method = result.data;
	    showPresets();
	},
	error => alert(reason)
    );
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
    for (let i=0; i<f.length; i++){
	status.innerHTML = " Loading file " + (i+1) + " of " + f.length;
	txt = await readFile(f[i]);
	shinylight.call('upload', {f: txt}, null).then(
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
