// Contains utility functions, classes, or helper methods that can be reused across your project.

import * as THREE from 'three';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import * as CANNON from 'cannon-es';

var threeScene;
var renderer;
var camera;
var cameraControl;
var grabber;
var mouseDown = false;
var objects = [];

const physicsScene = new CANNON.World({
    gravity: new CANNON.Vec3(0, -9.81, 0),
});

const timeStep = 1.0 / 60.0;

class Cloth {
    constructor(Nx, Ny, mass, clothSize){
        this.Nx = Nx;
        this.Ny= Ny;
        this.mass = mass;
        this.clothSize = clothSize;
        this.dist = clothSize / Nx;

        this.particles = [];
        this.shape = new CANNON.Particle();

        this.createParticles();
        this.connectParticles();
        this.createClothMesh();
    }

    createParticles(){
        const yOffset = 1.0;

        for(let i = 0; i < this.Nx + 1; i++) {
            this.particles.push([]);
            for(let j = 0; j < this.Ny + 1; j++) {
                const particle = new CANNON.Body({
                    mass: j === this.Ny ? 0 : this.mass,
                    shape: this.shape,
                    position: new CANNON.Vec3((i - this.Nx * 0.5) * this.dist, (j - this.Ny * 0.5) * this.dist + yOffset , 0),
                    velocity: new CANNON.Vec3(0, 0, -0.1 * (this.Ny - j))
                });
                this.particles[i].push(particle);
                physicsScene.addBody(particle);
            }
        }
    }

    connectParticles() {
        for(let i = 0; i < this.Nx + 1; i++) {
            for(let j = 0; j < this.Ny + 1; j++) {
                if(i < this.Nx)
                    this.connect(i, j, i + 1, j);
                if(j < this.Ny)
                    this.connect(i, j, i, j + 1);
            }
        }
    }

    connect(i1, j1, i2, j2) {
        physicsScene.addConstraint(new CANNON.DistanceConstraint(
            this.particles[i1][j1],
            this.particles[i2][j2],
            this.dist
        ));
    }

    createClothMesh() {
        this.clothGeometry = new THREE.PlaneGeometry(this.clothSize, this.clothSize, this.Nx, this.Ny);
        this.clothMat = new THREE.MeshPhongMaterial({
            side: THREE.DoubleSide,
            //wireframe: true,
            map: new THREE.TextureLoader().load('../assets/textures/adieu-ammenotep-1960.jpg')
        });

        this.clothMesh = new THREE.Mesh(this.clothGeometry, this.clothMat);
        threeScene.add(this.clothMesh);
    }

    simulate() {
        
        for(let i = 0; i < this.Nx + 1; i++) {
            for(let j = 0; j < this.Ny + 1; j++) {
                const index = j * (this.Nx + 1) + i;
    
                const positionAttribute = this.clothGeometry.attributes.position;
    
                const position = this.particles[i][this.Ny - j].position;
    
                positionAttribute.setXYZ(index, position.x, position.y, position.z);
    
                positionAttribute.needsUpdate = true;
            }
        }
        physicsScene.step(timeStep);
        renderer.render(threeScene, camera)

    }
}

class Ball {

    constructor(pos, radius, vel, world) {
        this.world = world;
        this.pos = pos;
        this.vel = vel;
        this.radius = radius;

        var geometry = new THREE.SphereGeometry( radius, 32, 32 )
        var material = new THREE.MeshPhongMaterial ({ color: 0xff0000});
        this.visMesh = new THREE.Mesh (geometry, material);
        threeScene.add(this.visMesh);
        this.visMesh.position.copy(pos);

        var sphereShape = new CANNON.Sphere(radius * 1.3);
        var sphereBody = new CANNON.Body({
            shape: sphereShape,
            position: pos,
            mass: 1
        });
        physicsScene.addBody(sphereBody);
    }

    simulate() {
        physicsScene.step(timeStep);
        renderer.render(threeScene, camera);

    }
}

export function initPhysScene() {

    var radius = 0.2;
    var pos = new THREE.Vector3(0, 1, 0);
    var vel = new THREE.Vector3(20.0, 20.0, 20.0);
    objects.push(new Ball(pos, radius, vel, physicsScene));

    var Nx = 15;
    var Ny = 15;
    var mass = 1;
    var clothSize = 1;
    objects.push(new Cloth(Nx, Ny, mass, clothSize))
}

export function initThreeScene() {
    threeScene = new THREE.Scene();
    
    threeScene.add( new THREE.AmbientLight( 0x505050 ) );	
    threeScene.fog = new THREE.Fog( 0x000000, 0, 15 );				

    var spotLight = new THREE.SpotLight( 0xffffff );
    spotLight.angle = Math.PI / 5;
    spotLight.penumbra = 0.1;
    spotLight.position.set( 0, 10, 0 );
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

    var ground = new THREE.Mesh(
        new THREE.PlaneGeometry( 20, 20, 1, 1 ),
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

    renderer = new THREE.WebGLRenderer();
    renderer.shadowMap.enabled = true;
    renderer.setPixelRatio( window.devicePixelRatio );
    renderer.setSize( 0.8 * window.innerWidth, 0.8 * window.innerHeight );

    window.addEventListener( 'resize', onWindowResize, false );
    container.appendChild( renderer.domElement );
            
    camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 0.01, 100);
    camera.position.set(0, 1, 4);
    camera.updateMatrixWorld();	
    threeScene.add( camera );

    cameraControl = new OrbitControls(camera, renderer.domElement);
    cameraControl.zoomSpeed = 2.0;
    cameraControl.panSpeed = 0.4;
}

window.addEventListener('resize', function() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
});

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize( window.innerWidth, window.innerHeight );
}

export function update() {
    for (var i = 0; i < objects.length; i++){
        objects[i].simulate();
    }
    renderer.render(threeScene, camera);
    cameraControl.update();

    requestAnimationFrame(update);
}

