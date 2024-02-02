import * as utils from './utils.js'

var resetPointer = utils.setupScene;
var runPointer = utils.run;
var stepPointer = utils.step;
var analyticSolutionPointer = utils.toggleSolution;

var resetButton = document.getElementById('resetButton');
var runButton = document.getElementById('runButton');
var stepButton = document.getElementById('stepButton');
var toggleAnalyticSolution = document.getElementById('toggleAnalyticSolution');

resetButton.onclick = resetPointer;
runButton.onclick = runPointer;
stepButton.onclick = stepPointer;
toggleAnalyticSolution.onclick = analyticSolutionPointer;

utils.setupScene();
utils.update();	
