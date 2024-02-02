import * as utils from "./utils.js"

var runPointer = utils.run;
var restartPointer = utils.restart;
var edgePointer = utils.onShowEdges;
var collisionPointer = utils.onCollision;

var runButton = document.getElementById("runButtonId");
var restartButton = document.getElementById("restartButtonId");
var edgeButton = document.getElementById("edgeButtonId");
var collisionButton = document.getElementById("collisionButtonId");

runButton.onclick = runPointer;
restartButton.onclick = restartPointer;
edgeButton.onclick = edgePointer;
collisionButton.onclick = collisionPointer;

utils.initThreeScene();
utils.onWindowResize();
utils.initPhysics();
utils.update();
			