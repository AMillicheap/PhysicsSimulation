import * as utils from "./utils.js"

var restartPointer = utils.setupScene;
var restartButton = document.getElementById("restartButtonId");
restartButton.onclick = restartPointer;

utils.setupScene();
utils.update();
