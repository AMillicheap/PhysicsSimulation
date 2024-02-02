var canvas = document.getElementById("myCanvas");
var c = canvas.getContext("2d");
canvas.width = window.innerWidth - 20;
canvas.height = window.innerHeight - 20;

var simMinWidth = 1.0;
var cScale = 0.014 * Math.min(canvas.width, canvas.height) / simMinWidth;

function cX(pos) { return canvas.width / 2 + pos.x * cScale; }
function cY(pos) { return 0.4 * canvas.height - pos.y * cScale; }

class Vector2 {
    constructor(x = 0.0, y = 0.0) {
        this.x = x; 
        this.y = y;
    }

    set(v) {
        this.x = v.x; this.y = v.y;
    }

    clone() {
        return new Vector2(this.x, this.y);
    }

    add(v, s = 1.0) {
        this.x += v.x * s;
        this.y += v.y * s;
        return this;
    }

    addVectors(a, b) {
        this.x = a.x + b.x;
        this.y = a.y + b.y;
        return this;
    }

    subtract(v, s = 1.0) {
        this.x -= v.x * s;
        this.y -= v.y * s;
        return this;
    }

    subtractVectors(a, b) {
        this.x = a.x - b.x;
        this.y = a.y - b.y;
        return this;			
    }

    length() {
        return Math.sqrt(this.x * this.x + this.y * this.y);
    }

    scale(s) {
        this.x *= s;
        this.y *= s;
        return this;
    }

    dot(v) {
        return this.x * v.x + this.y * v.y;
    }

    perp() {
        return new Vector2(-this.y, this.x);
    }
}

class Particle {
    constructor(pos, mass) {
        this.pos = pos.clone();
        this.prevPos = pos.clone();
        this.vel = new Vector2();
        this.mass = mass;
    }

    startStep(dt, gravity) {
        this.vel.add(gravity, dt);
        this.prevPos.set(this.pos);
        this.pos.add(this.vel, dt);
    }

    endStep(dt) {
        this.vel.subtractVectors(this.pos, this.prevPos);
        this.vel.scale(1.0 / dt);
    }
}
class Hair {
    constructor(color, angles) {
        this.color = color;
        this.masses = [0.0];
        this.pos = [{x:0.0, y:0.0}];
        this.prevPos = [{x:0.0, y:0.0}];
        this.vel = [{x:0.0, y:0.0}];
        this.theta = [0.0];
        this.omega = [0.0];
        this.particleDist = 0.06

        this.trail = new Int32Array(1000);
        this.trailFirst = 0;
        this.trailLast = 0;

        var x = 0.0, y = 0.0;
        for (var i = 0; i < angles.length; i++) {
            this.masses.push(1.0);
            this.theta.push(angles[i]);
            this.omega.push(0.0);
            x += this.particleDist * Math.sin(angles[i]);
            y += this.particleDist * -Math.cos(angles[i]); 
            this.pos.push({ x:x, y:y});
            this.prevPos.push({ x:x, y:y});
            this.vel.push({x:0, y:0});
        }
    }

    simulate(dt, gravity, dampingFactor) 
    {
        var p = this;
        for (var i = 1; i < p.masses.length; i++) {
            p.vel[i].y += dt * scene.gravity;
            p.prevPos[i].x = p.pos[i].x;
            p.prevPos[i].y = p.pos[i].y;
            p.pos[i].x += p.vel[i].x * dt;
            p.pos[i].y += p.vel[i].y * dt;
        }
        for (var i = 1; i < p.masses.length; i++) {
            var dx = p.pos[i].x - p.pos[i-1].x;
            var dy = p.pos[i].y - p.pos[i-1].y;
            var d = Math.sqrt(dx * dx + dy * dy);
            var w0 = p.masses[i - 1] > 0.0 ? 1.0 / p.masses[i - 1] : 0.0;
            var w1 = p.masses[i] > 0.0 ? 1.0 / p.masses[i] : 0.0;
            var corr = (p.particleDist - d) / d / (w0 + w1);
            p.pos[i - 1].x -= w0 * corr * dx; 
            p.pos[i - 1].y -= w0 * corr * dy; 
            p.pos[i].x += w1 * corr * dx; 
            p.pos[i].y += w1 * corr * dy; 
        }
        for (var i = 1; i < p.masses.length; i++) {
            p.vel[i].x = (p.pos[i].x - p.prevPos[i].x) / dt;
            p.vel[i].y = (p.pos[i].y - p.prevPos[i].y) / dt;
            var dampingForce = { x: -dampingFactor * p.vel[i].x, y: -dampingFactor * p.vel[i].y };
            p.vel[i].x += dampingForce.x * dt;
            p.vel[i].y += dampingForce.y * dt;
        }
    }
    
    updateTrail() {
        this.trail[this.trailLast] = cX(this.pos[this.pos.length-1]);
        this.trail[this.trailLast + 1] = cY(this.pos[this.pos.length-1]);
        this.trailLast = (this.trailLast + 2) % this.trail.length;
        if (this.trailLast == this.trailFirst)
            this.trailFirst = (this.trailFirst + 2) % this.trail.length;
    }
    draw() {
        c.strokeStyle = this.color;
        c.lineWidth = 2.0;
        if (this.trailLast != this.trailFirst) {
            var i = this.trailFirst;
            c.beginPath();
            c.moveTo(this.trail[i], this.trail[i + 1]);
            i = (i + 2) % this.trail.length;
            while (i != this.trailLast) {
                c.lineTo(this.trail[i], this.trail[i + 1]);
                i = (i + 2) % this.trail.length;
            }
            c.stroke();
        }

        var p = this;
        c.strokeStyle = "#303030";
        c.lineWidth = 10;
        c.beginPath();
        c.moveTo(cX(p.pos[0]), cY(p.pos[0]));
        console.log(p);
        for (var i = 1; i < p.masses.length; i++) 
            c.lineTo(cX(p.pos[i]), cY(p.pos[i]));
        c.stroke();
        c.lineWidth = 1;            

        c.fillStyle = this.color;
        for (var i = 1; i < p.masses.length; i++) {
            var r = 0.03 * Math.sqrt(p.theta[i]);
            c.beginPath();			
            c.arc(
                cX(p.pos[i]), cY(p.pos[i]), 2 * cScale * r, 0.0, 2.0 * Math.PI); 
            c.closePath();
            c.fill();
        }
    }
}

var scene = {
    gravity : -10.0,
    dt : 0.01,
    numSubSteps : 10000,
    paused : true,
    pendulumPBD : null,
    dampingFactor : 0.15
};

var sceneNr = 0;

export function setupScene() {

    var numParticles = 500;
    var angles = [0.5 * Math.PI];
    var masses = [0];

    for (var i = 1; i < numParticles; i++) {
        angles.push(0.5 * Math.PI);
        masses.push(1.0);
    }

    scene.pendulumPBD = new Hair("#FF3030", angles);

    scene.paused = true;
    sceneNr++;
}

function draw() {
    c.fillStyle = "#000000";
    c.fillRect(0, 0, canvas.width, canvas.height);
    scene.pendulumPBD.draw();
}

function simulate() {
    if (scene.paused)
        return;
    var sdt = scene.dt / scene.numSubSteps;
    var trailGap = scene.numSubSteps / 10;

    for (var step = 0; step < scene.numSubSteps; step++) {
        scene.pendulumPBD.simulate(sdt, scene.gravity, scene.dampingFactor);
        if (step % trailGap == 0)
            scene.pendulumPBD.updateTrail();
    }
}

document.getElementById("stepsSlider").oninput = function() {
    var steps = [1, 5, 10, 100, 1000, 10000];
    scene.numSubSteps = steps[Number(this.value)];
    document.getElementById("steps").innerHTML = scene.numSubSteps.toString();
}

document.addEventListener("keydown", event => {
    if (event.isComposing || event.keyCode === 229) 
        return;
    if (event.key == 's')
        step();
    });    

export function run() {
    scene.paused = false;
}

function step() {
    scene.paused = false;
    simulate();
    scene.paused = true;
}

export function update() {
    simulate();
    draw();
    requestAnimationFrame(update);
}