import * as utils from "./utils.js"

var runPointer = utils.run;
var restartPointer = utils.restart;
var squashPointer = utils.squash;
var newBodyPointer = utils.newBody;

var runButton = document.getElementById("buttonRun");
var restartButton = document.getElementById("restartButton");
var squashButton = document.getElementById("squashButton");
var newBodyButton = document.getElementById("newBodyButton")

runButton.onclick = runPointer;
restartButton.onclick = restartPointer;
squashButton.onclick = squashPointer;
newBodyButton.onclick = newBodyPointer;

utils.initThreeScene();
utils.onWindowResize();
utils.initPhysics();
utils.update();
