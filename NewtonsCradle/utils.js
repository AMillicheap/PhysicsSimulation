// drawing -------------------------------------------------------

var canvas = document.getElementById("myCanvas");
var c = canvas.getContext("2d");

canvas.width = window.innerWidth - 20;
canvas.height = window.innerHeight - 100;

var simMinWidth = 2.0;
var cScale = Math.min(canvas.width, canvas.height) / simMinWidth;
var simWidth = canvas.width / cScale;
var simHeight = canvas.height / cScale;

function cX(pos) {
	return pos.x * cScale;
}

function cY(pos) {
	return canvas.height - pos.y * cScale;
}

// vector math -------------------------------------------------------

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

// scene -------------------------------------------------------

class Bead {
	constructor(radius, mass, pos) {
		this.radius = radius;
		this.mass = mass;
		this.pos = pos.clone();
		this.prevPos = pos.clone();
		this.vel = new Vector2();
	}
	startStep(dt, gravity) {
		this.vel.add(gravity, dt);
		this.prevPos.set(this.pos);
		this.pos.add(this.vel, dt);
	}
	keepOnWire(center, radius) {
		var dir = new Vector2();
		dir.subtractVectors(this.pos, center);
		var len = dir.length();
		if (len == 0.0)
			return;
		dir.scale(1.0 / len);
		var lambda = physicsScene.wireRadius - len;
		this.pos.add(dir, lambda);
		return lambda;
	}
	endStep(dt) {
		this.vel.subtractVectors(this.pos, this.prevPos);
		this.vel.scale(1.0 / dt);
	}

}

var physicsScene = 
{
	gravity : new Vector2(0.0, -10.0),
	dt : 1.0 / 60.0,
	numSteps : 100,
	wireCenter : new Vector2(),
	wireRadius : 0.0,
	ballRadius : 0.1,
	beads : [],
	wires : [],
	selectedBead : null,
	dampingFactor : 0.5,
};

// -----------------------------------------------------

export function setupScene() 
{
	physicsScene.beads = [];
	physicsScene.wires = [];
	var numBeads = 5;

	physicsScene.wireRadius = simMinWidth * 0.3;
	var mass = 1.0;
	var angle = 3 * Math.PI / 2;

	for (var i = 0; i < numBeads; i++) {
		var x = (simWidth / 2.0) - (2 * (2 - i) * physicsScene.ballRadius);
		var y = (simHeight / 2.0);
		physicsScene.wires.push(new Vector2(x, y));
	}

	for (var i = 0; i < numBeads; i++) {
		var mass = Math.PI * physicsScene.ballRadius * physicsScene.ballRadius;			
		var pos = new Vector2(
			physicsScene.wires[i].x + physicsScene.wireRadius * Math.cos(angle), 
			physicsScene.wires[i].y + physicsScene.wireRadius * Math.sin(angle));

		physicsScene.beads.push(new Bead(physicsScene.ballRadius, mass, pos));
	}
}

// draw -------------------------------------------------------

function drawCircle(pos, radius, filled)
{
	c.beginPath();			
	c.arc(
		cX(pos), cY(pos), cScale * radius, 0.0, 2.0 * Math.PI); 
	c.closePath();
	if (filled)
		c.fill();
	else 
		c.stroke();
}

function drawRect(pos, width, height, filled) {
	var pos1 = new Vector2(pos.x + width, pos.y);
	var pos2 = new Vector2(pos.x + width, pos.y + height);
	var pos3 = new Vector2(pos.x, pos.y + height);
	c.beginPath();
	c.moveTo(cX(pos), cY(pos));
	c.lineTo(cX(pos1), cY(pos1));
	c.lineTo(cX(pos2), cY(pos2));
	c.lineTo(cX(pos3), cY(pos3));
	c.closePath();
	c.stroke();
	if (filled) {
		c.fill();
	}
}

function draw() 
{
	c.clearRect(0, 0, canvas.width, canvas.height);

	c.fillStyle = "#000000";
	c.lineWidth = 2.0;

	
	var posRect = new Vector2(physicsScene.wires[0].x - 0.1 * simMinWidth, physicsScene.wires[0].y);
	drawRect(posRect, 9 * physicsScene.ballRadius + 0.2 * simMinWidth, 0.02 * simMinWidth)
	
	c.font = "48px serif";
	c.fillText("Move the ball bearings manually", 420, 50);

	c.fillStyle = "#FF0000"
	for (var i = 0; i < physicsScene.beads.length; i++) {
		var bead = physicsScene.beads[i];
		drawCircle(bead.pos, bead.radius, true);
		c.beginPath();
		c.moveTo(cX(physicsScene.wires[i]), cY(physicsScene.wires[i]));
		c.lineTo(cX(bead.pos), cY(bead.pos));
		c.stroke();
	}
}

// --- collision handling -------------------------------------------------------

function handleBeadBeadCollision(bead1, bead2) 
{
	var restitution = 1.0;
	var dir = new Vector2();
	dir.subtractVectors(bead2.pos, bead1.pos);
	var d = dir.length();
	if (d == 0.0 || d > bead1.radius + bead2.radius)
		return;

	dir.scale(1.0 / d);

	var corr = (bead1.radius + bead2.radius - d) / 2.0;
	bead1.pos.add(dir, -corr);
	bead2.pos.add(dir, corr);

	var v1 = bead1.vel.dot(dir);
	var v2 = bead2.vel.dot(dir);

	var m1 = bead1.mass;
	var m2 = bead2.mass;

	var newV1 = (m1 * v1 + m2 * v2 - m2 * (v1 - v2) * restitution) / (m1 + m2);
	var newV2 = (m1 * v1 + m2 * v2 - m1 * (v2 - v1) * restitution) / (m1 + m2);

	bead1.vel.add(dir, newV1 - v1);
	bead2.vel.add(dir, newV2 - v2);
}

canvas.addEventListener('mousedown', function (event) {
	var mouseX = event.clientX - canvas.getBoundingClientRect().left;
	var mouseY = event.clientY - canvas.getBoundingClientRect().top;

	for (var i = 0; i < physicsScene.beads.length; i++) {
		var bead = physicsScene.beads[i];
		var distance = Math.sqrt(Math.pow(cX(bead.pos) - mouseX, 2) + Math.pow(cY(bead.pos) - mouseY, 2));

		if (distance <= cScale * bead.radius) {
			physicsScene.selectedBead = bead;
			break;
		}
	}
});

canvas.addEventListener('mousemove', function (event) {
	if (physicsScene.selectedBead) {
		var mouseX = event.clientX - canvas.getBoundingClientRect().left;
		var mouseY = event.clientY - canvas.getBoundingClientRect().top;

		physicsScene.selectedBead.pos.x = mouseX / cScale;
		physicsScene.selectedBead.pos.y = (canvas.height - mouseY) / cScale;
		physicsScene.selectedBead.vel.scale(physicsScene.dampingFactor);
	}
});

canvas.addEventListener('mouseup', function (event) {
	physicsScene.selectedBead = null;
});


function simulate() 
{
	var sdt = physicsScene.dt / physicsScene.numSteps;

	for (var step = 0; step < physicsScene.numSteps; step++) {
		for (var i = 0; i < physicsScene.beads.length; i++)
			physicsScene.beads[i].startStep(sdt, physicsScene.gravity);

		for (var i = 0; i < physicsScene.beads.length; i++) {
			physicsScene.beads[i].keepOnWire(
				physicsScene.wires[i], physicsScene.wireRadius);
		}

		for (var i = 0; i < physicsScene.beads.length; i++)
			physicsScene.beads[i].endStep(sdt);

		for (var i = 0; i < physicsScene.beads.length; i++) {
			for (var j = 0; j < i; j++) {
				handleBeadBeadCollision(
					physicsScene.beads[i], physicsScene.beads[j]);
			}
		}
	}
}

// --------------------------------------------------------

export function update() {
	simulate();
	draw();
	requestAnimationFrame(update);
}