import * as utils from './utils.js' 

var runPointer = utils.run;
var restartPointer = utils.restart;
var collisionPointer = utils.onShowColl;

var runButton = document.getElementById("buttonRun");
var restartButton = document.getElementById("buttonRestart");
var collisionButton = document.getElementById("checkColl");

runButton.onclick = runPointer;
restartButton.onclick = restartPointer;
collisionButton.onclick = collisionPointer;

utils.initThreeScene();
utils.onWindowResize();
utils.initPhysics();
utils.update();
