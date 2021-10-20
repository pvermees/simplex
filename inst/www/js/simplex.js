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

function loaddefaults(fn){
    shinylight.call('loadTable', {fn: fn}, null).then(function(result) {
	result.data = result.data.tab;
        shinylight.setGridResult(inputTable, result);
    }).catch(function(reason) {
        shinylight.setElementText('error', reason);
    });
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

function setup(){
    selectButton(0);
    loadPage("setup.html")
	.then(
	    () => { if (simplex.method==null) { alert('empty') } }
	);
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
