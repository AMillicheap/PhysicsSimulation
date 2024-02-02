import * as utils from './utils.js'

var runPointer = utils.run;
var restartPointer = utils.setupScene;
var stepPointer = utils.step;

var runButton = document.getElementById("runButton");
var resetButton = document.getElementById("resetButton");
var stepButton = document.getElementById("stepButton");

runButton.onclick = runPointer;
resetButton.onclick = restartPointer;
stepButton.onclick = stepPointer;

utils.setupScene();
utils.update();	
