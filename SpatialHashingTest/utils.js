
var threeScene;
var renderer;
var camera;
var cameraControl;	

// ------------------------------------------------------------------
class Hash {
    constructor(spacing, maxNumObjects) 
    {
        this.spacing = spacing;
        this.tableSize = 2 * maxNumObjects;

        // cellStart stores the starting index for each hash cell in the the hash table.
        // cellEntries stores the indices of objects within the hash cells.

        this.cellStart = new Int32Array(this.tableSize + 1);
        this.cellEntries = new Int32Array(maxNumObjects);
        this.queryIds = new Int32Array(maxNumObjects);
        this.querySize = 0;
    }

    hashCoords(xi, yi, zi) {
        var h = (xi * 92837111) ^ (yi * 689287499) ^ (zi * 283923481);	// fantasy function
        return Math.abs(h) % this.tableSize; 
    }

    intCoord(coord) {
        return Math.floor(coord / this.spacing);
    }

    hashPos(pos, nr) {
        return this.hashCoords(
            this.intCoord(pos[3 * nr]), 
            this.intCoord(pos[3 * nr + 1]),
            this.intCoord(pos[3 * nr + 2]));
    }

    create(pos) {
        var numObjects = Math.min(this.maxNumObjects, this.cellEntries.length);

        this.cellStart.fill(0);
        this.cellEntries.fill(0);

        for (var i = 0; i < numObjects; i++) {
            var h = this.hashPos(pos, i);
            this.cellStart[h]++;
        }

        var start = 0;
        for (var i = 0; i < this.tableSize; i++) {
            start += this.cellStart[i];
            this.cellStart[i] = start;
        }
        this.cellStart[this.tableSize] = start;

        for (var i = 0; i < numObjects; i++) {
            var h = this.hashPos(pos, i);
            this.cellStart[h]--;
            this.cellEntries[this.cellStart[h]] = i;
        }
    }

    query(pos, nr, maxDist) {
        var x0 = this.intCoord(pos[3 * nr] - maxDist);
        var y0 = this.intCoord(pos[3 * nr + 1] - maxDist);
        var z0 = this.intCoord(pos[3 * nr + 2] - maxDist);

        var x1 = this.intCoord(pos[3 * nr] + maxDist);
        var y1 = this.intCoord(pos[3 * nr + 1] + maxDist);
        var z1 = this.intCoord(pos[3 * nr + 2] + maxDist);

        this.querySize = 0;

        for (var xi = x0; xi <= x1; xi++) {
            for (var yi = y0; yi <= y1; yi++) {
                for (var zi = z0; zi <= z1; zi++) {
                    var h = this.hashCoords(xi, yi, zi);
                    var start = this.cellStart[h];
                    var end = this.cellStart[h + 1];

                    for (var i = start; i < end; i++) {
                        this.queryIds[this.querySize] = this.cellEntries[i];
                        this.querySize++;
                    }
                }
            }
        }
    }
};

// ----- math on vector arrays -------------------------------------------------------------

function vecScale(a,anr, scale) {
    anr *= 3;
    a[anr++] *= scale;
    a[anr++] *= scale;
    a[anr]   *= scale;
}

function vecCopy(a,anr, b,bnr) {
    anr *= 3; bnr *= 3;
    a[anr++] = b[bnr++]; 
    a[anr++] = b[bnr++]; 
    a[anr]   = b[bnr];
}

function vecAdd(a,anr, b,bnr, scale = 1.0) {
    anr *= 3; bnr *= 3;
    a[anr++] += b[bnr++] * scale; 
    a[anr++] += b[bnr++] * scale; 
    a[anr]   += b[bnr] * scale;
}

function vecSetDiff(dst,dnr, a,anr, b,bnr, scale = 1.0) {
    dnr *= 3; anr *= 3; bnr *= 3;
    dst[dnr++] = (a[anr++] - b[bnr++]) * scale;
    dst[dnr++] = (a[anr++] - b[bnr++]) * scale;
    dst[dnr]   = (a[anr] - b[bnr]) * scale;
}

function vecLengthSquared(a,anr) {
    anr *= 3;
    let a0 = a[anr], a1 = a[anr + 1], a2 = a[anr + 2];
    return a0 * a0 + a1 * a1 + a2 * a2;
}

function vecDistSquared(a,anr, b,bnr) {
    anr *= 3; bnr *= 3;
    let a0 = a[anr] - b[bnr], a1 = a[anr + 1] - b[bnr + 1], a2 = a[anr + 2] - b[bnr + 2];
    return a0 * a0 + a1 * a1 + a2 * a2;
}	

function vecDot(a,anr, b,bnr) {
    anr *= 3; bnr *= 3;
    return a[anr] * b[bnr] + a[anr + 1] * b[bnr + 1] + a[anr + 2] * b[bnr + 2];
}	

// ------------------------------------------------------------------

var physicsScene = 
{
    gravity : [0.0, 0.0, 0.0],
    dt : 1.0 / 60.0,
    worldBounds :  [-1.0, 0.0, -1.0, 1.0, 2.0, 1.0],
    paused: true,
    balls: null,
};

export function onShowColl() {
    if (physicsScene.balls)
        physicsScene.balls.showCollisions = !physicsScene.balls.showCollisions;
}			
            
// ------------------------------------------------------------------
class Balls {
    constructor(radius, pos, vel, numBalls, scene)
    {
        this.radius = radius;
        this.pos = pos;
        this.prevPos = pos;
        this.vel = vel;
        this.matrix = new THREE.Matrix4();
        this.numBalls = numBalls;
        this.hash = new Hash(2.0 * radius, this.numBalls);
        this.showCollisions = false;

        this.normal = new Float32Array(3);

        // visual mesh

        var geometry = new THREE.SphereGeometry( radius, 8, 8 );
        var material = new THREE.MeshPhongMaterial();

        this.visMesh = new THREE.InstancedMesh( geometry, material, this.numBalls );
        this.visMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage); 

        this.ballColor = new THREE.Color(0x000000);
        this.ballCollisionColor = new THREE.Color(0x39FF14);

        var colors = new Float32Array(3 * this.numBalls);
        this.visMesh.instanceColor = new THREE.InstancedBufferAttribute(colors, 3, false, 1);
        for (var i = 0; i < this.numBalls; i++) 
            this.visMesh.setColorAt(i, this.ballColor);

        threeScene.add(this.visMesh);

        this.updateMesh();
    }

    updateMesh()
    {
        for (var i = 0; i < this.numBalls; i++) {
            this.matrix.makeTranslation(this.pos[3 * i], this.pos[3 * i + 1], this.pos[3 * i + 2]);
            this.visMesh.setMatrixAt(i, this.matrix);
        }
        this.visMesh.instanceMatrix.needsUpdate = true;
        this.visMesh.instanceColor.needsUpdate = true;
    }

    simulate(dt, gravity, worldBounds)
    {
        var minDist = 2.0 * this.radius;

        // integrate

        for (var i = 0; i < this.numBalls; i++) {
            vecAdd(this.vel, i, gravity, 0, dt);
            vecCopy(this.prevPos, i, this.pos, i);
            vecAdd(this.pos, i, this.vel, i, dt);
        }

        this.hash.create(this.pos);

        // handle collisions

        for (var i = 0; i < this.numBalls; i++) {

            this.visMesh.setColorAt(i, this.ballColor);

            // world collision

            for (var dim = 0; dim < 3; dim++) {

                var nr = 3 * i + dim;
                if (this.pos[nr] < worldBounds[dim] + this.radius) {
                    this.pos[nr] = worldBounds[dim] + this.radius;
                    this.vel[nr] = - this.vel[nr];
                    if (this.showCollisions)
                        this.visMesh.setColorAt(i, this.ballCollisionColor);
                }
                else if (this.pos[nr] > worldBounds[dim + 3] - this.radius) {
                    this.pos[nr] = worldBounds[dim + 3] - this.radius;
                    this.vel[nr] = - this.vel[nr];
                    if (this.showCollisions)
                        this.visMesh.setColorAt(i, this.ballCollisionColor);
                }
            }

            //  interball collision

            this.hash.query(this.pos, i, 2.0 * this.radius);

            for (var nr = 0; nr < this.hash.querySize; nr++) {
                var j = this.hash.queryIds[nr];

                vecSetDiff(this.normal, 0, this.pos, i, this.pos, j);
                var d2 = vecLengthSquared(this.normal, 0);

                 // are the balls overlapping?

                if (d2 > 0.0 && d2 < minDist * minDist) {
                    var d = Math.sqrt(d2);
                    vecScale(this.normal, 0, 1.0 / d);	

                    // separate the balls

                    var corr = (minDist - d) * 0.5;

                    vecAdd(this.pos, i, this.normal, 0, corr);
                    vecAdd(this.pos, j, this.normal, 0, -corr);

                    // reflect velocities along normal

                    var vi = vecDot(this.vel, i, this.normal, 0);
                    var vj = vecDot(this.vel, j, this.normal, 0);

                    vecAdd(this.vel, i, this.normal, 0, vj - vi);
                    vecAdd(this.vel, j, this.normal, 0, vi - vj);

                    if (this.showCollisions)
                        this.visMesh.setColorAt(i, this.ballCollisionColor);
                }
            }
        }
        this.updateMesh();
    }
}

var initialNumBalls = parseInt(document.getElementById("particleCountSlider").value);

// ------------------------------------------------------------------
export function initPhysics(scene) 
{
    var radius = 0.025;

    var spacing = 3.0 * radius;
    var velRand = 0.2;

    var s = physicsScene.worldBounds;

    var numX = Math.floor((s[3] - s[0] - 2.0 * spacing) / spacing);
    var numY = Math.floor((s[4] - s[1] - 2.0 * spacing) / spacing);
    var numZ = Math.floor((s[5] - s[2] - 2.0 * spacing) / spacing);

    var pos = new Float32Array(3 * numX * numY * numZ);
    var vel = new Float32Array(3 * numX * numY * numZ);
    vel.fill(0.0);

    for (var xi = 0; xi < numX; xi++) {
        for (var yi = 0; yi < numY; yi++) {
            for (var zi = 0; zi < numZ; zi++) {
                var x = 3 * ((xi * numY + yi) * numZ + zi);
                var y = x + 1;
                var z = x + 2;
                pos[x] = s[0] + spacing + xi * spacing;
                pos[y] = s[1] + spacing + yi * spacing;
                pos[z] = s[2] + spacing + zi * spacing;

                vel[x] = -velRand + 2.0 * velRand * Math.random();
                vel[y] = -velRand + 2.0 * velRand * Math.random();
                vel[z] = -velRand + 2.0 * velRand * Math.random();
            }
        }
    }

    physicsScene.balls = new Balls(radius, pos, vel, initialNumBalls, threeScene);

    //document.getElementById("particleCount").innerHTML = pos.length / 3;		

}

document.getElementById("particleCountSlider").addEventListener("input", function() {
    // Update numBalls based on the slider value
    physicsScene.balls.numBalls = parseInt(this.value);
    // Update the label to display the current value
    document.getElementById("particleCountLabel").textContent = this.value;
    // Update the Balls instance
    physicsScene.balls.updateMesh();
});

var timeFrames = 0;
var timeSum = 0;	
            
// ------------------------------------------------------------------
function simulate() 
{
    if (physicsScene.paused)
        return;

    var startTime = performance.now();					

    physicsScene.balls.simulate(physicsScene.dt, 
        physicsScene.gravity, physicsScene.worldBounds);

    var endTime = performance.now();
    timeSum += endTime - startTime; 
    timeFrames++;

    if (timeFrames > 10) {
        timeSum /= timeFrames;
        document.getElementById("ms").innerHTML = timeSum.toFixed(3);		
        timeFrames = 0;
        timeSum = 0;
    }
}
            
// ------------------------------------------
        
export function initThreeScene() 
{
    threeScene = new THREE.Scene();
    
    // Lights
    
    threeScene.add( new THREE.AmbientLight( 0x505050 ) );	
    threeScene.fog = new THREE.Fog( 0x000000, 0, 15 );				

    var spotLight = new THREE.SpotLight( 0xffffff );
    spotLight.angle = Math.PI / 5;
    spotLight.penumbra = 0.2;
    spotLight.position.set( 2, 3, 3 );
    spotLight.castShadow = true;
    spotLight.shadow.camera.near = 3;
    spotLight.shadow.camera.far = 10;
    spotLight.shadow.mapSize.width = 1024;
    spotLight.shadow.mapSize.height = 1024;
    threeScene.add( spotLight );

    var dirLight = new THREE.DirectionalLight( 0x55505a, 1 );
    dirLight.position.set( 0, 3, 0 );
    dirLight.castShadow = true;
    dirLight.shadow.camera.near = 1;
    dirLight.shadow.camera.far = 10;

    dirLight.shadow.camera.right = 1;
    dirLight.shadow.camera.left = - 1;
    dirLight.shadow.camera.top	= 1;
    dirLight.shadow.camera.bottom = - 1;

    dirLight.shadow.mapSize.width = 1024;
    dirLight.shadow.mapSize.height = 1024;
    threeScene.add( dirLight );
    
    // Geometry

    var ground = new THREE.Mesh(
        new THREE.PlaneBufferGeometry( 20, 20, 1, 1 ),
        new THREE.MeshPhongMaterial( { color: 0xa0adaf, shininess: 150 } )
    );				

    ground.rotation.x = - Math.PI / 2; // rotates X/Y to X/Z
    ground.receiveShadow = true;
    threeScene.add( ground );
    
    var helper = new THREE.GridHelper( 20, 20 );
    helper.material.opacity = 1.0;
    helper.material.transparent = true;
    helper.position.set(0, 0.002, 0);
    threeScene.add( helper );				
    
    // Renderer

    renderer = new THREE.WebGLRenderer();
    renderer.shadowMap.enabled = true;
    renderer.setPixelRatio( window.devicePixelRatio );
    renderer.setSize( 0.8 * window.innerWidth, 0.8 * window.innerHeight );
    window.addEventListener( 'resize', onWindowResize, false );
    container.appendChild( renderer.domElement );
    
    // Camera
            
    camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 0.01, 100);
    camera.position.set(0, 2, 4);
    camera.updateMatrixWorld();	

    threeScene.add( camera );

    cameraControl = new THREE.OrbitControls(camera, renderer.domElement);
    cameraControl.zoomSpeed = 2.0;
    cameraControl.panSpeed = 0.4;
}

export function onWindowResize() {

    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize( window.innerWidth, window.innerHeight );
}

export function run() {
    var button = document.getElementById('buttonRun');
    if (physicsScene.paused)
        button.innerHTML = "Stop";
    else
        button.innerHTML = "Run";
    physicsScene.paused = !physicsScene.paused;
}

export function restart() {
    location.reload();
}

// make browser to call us repeatedly -----------------------------------

export function update() {
    simulate();
    renderer.render(threeScene, camera);
    cameraControl.update();				
    
    requestAnimationFrame(update);
}
