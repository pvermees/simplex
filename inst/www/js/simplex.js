var simplex = {
    'selected': 0,
    'buttonIDs': ['setup','drift','logratios','calibration','samples','finish'],
    'oncol': 'yellow',
    'offcol': 'white',
    'method': null
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

async function loadPage(url) {
    fetch(url)
	.then(response => response.text())
	.then(text => {
	    document.getElementById("contents").innerHTML = text;
	});
}

async function setup(){
    selectButton(0);
    loadPage("setup.html"); // TODO: load presets when page has loaded
}

function loadPresets(){
    let m = document.getElementById("methods").value;
    shinylight.call('presets', {method: m}, null).then(function(result) {
	simplex.method = result.data;
	showPresets();
    }).catch(function(reason) {
        alert(reason);
    });
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

function chooseFile(t){
    const selectedFile = document.getElementById(t.id).files[0];
    alert(selectedFile.length)
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
