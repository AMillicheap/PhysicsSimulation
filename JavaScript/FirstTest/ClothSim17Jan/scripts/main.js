// Initialize your Three.js scene, set up event listeners, and handle overall logic.

import * as utils from './utils.js';


utils.initThreeScene();
utils.initPhysScene();
utils.update();

window.addEventListener('resize', function() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
});