#include <SFML\Main.hpp>
#include <SFML\Graphics.hpp>
#include <vector>
#include <iostream>

//SETTINGS
bool GRAVITY_VERT_ENABLE = true;
bool GRAVITY_BETWEEN_BALLS = true;
const double BALL_RAD_START = 10;
const double BALL_RAD_MAX = 20;
const double BALL_RAD_MIN = 5;


//CONSTANTS
const double G = 0.67;
const double g = 9.807;
const double PI = 3.1415;
const double dt = 0.1;
const double collision_extra_dist = 0.08;
const double ball_density = 2;
const double ball_bouncyness = 0.6;
const double wall_bouncyness = 0.8;
const double slingShotStrength = 0.4;
const double PRED_VECTORS_START_ALPHA = 200;
const int PRED_VECTORS_AMOUNT = 100;
const double PRED_VECTOR_TIME_STEP = 0.1;

const double GRAVITY_WELL_CUTOFF_CLOSE_DIST = 15;
const double GRAVITY_FIELD_LINE_START_NUMBER = 20;



//TOGGLES & VALUES
sf::Vector3f newObjectInfoVector(0, 0, 0);
sf::Vector2f lastMouseFreeHand;
double newBallRad = BALL_RAD_START;
bool PAUSE = false;

//ELEMENTS OF THE SIM
struct Ball
{

	sf::CircleShape body;
	sf::Vector2f velocity;
	sf::Vector2f acceleration;

	Ball(sf::Vector2f p, double r)
	{
		acceleration = sf::Vector2f(0, 0);
		velocity = sf::Vector2f(0, 0);

		body.setRadius(r);
		body.setOrigin(r, r);
		body.setOutlineThickness(-1);
		body.setFillColor(sf::Color(0, 0, 0, 0));
		body.setOutlineColor(sf::Color::White);
		body.setPointCount(200);
		body.setPosition(p);

	}
	double mass() { return ball_density * body.getRadius() * body.getRadius(); }
	
};
std::vector<Ball> ballList;

struct Wall
{
	
	sf::Vector2f pos;
	double angle;
	double length;
	double killer;

	Wall(sf::Vector2f p,double l, double a, bool k)
	{
		pos = p;
		angle = a;
		length = l;
		killer = k;
	}

};
std::vector<Wall> wallList;

struct ballSpawner
{
	sf::Vector2f pos;
	sf::Vector2f spawnVel;
	double sizeBallSpawn;
	double delay;
	double timeSinceLastSpawn;

	ballSpawner(sf::Vector2f p, sf::Vector2f v, double size, double d)
	{
		pos = p, spawnVel = v, delay = d, timeSinceLastSpawn = 0, sizeBallSpawn = size;
	}

};
std::vector<ballSpawner> spawnerList;

struct gravityWell
{
	sf::Vector2f pos;
	double strength;
	double rad;

	gravityWell(sf::Vector2f p, double s, double r) { pos = p; strength = s; rad = r; }

};
std::vector<gravityWell> gravityWellList;


//||||||||||||FUNCTIONS||||||||||||||||||||


/*BALLS*/
void gravityVerticalBalls(std::vector<Ball> &l)
{
	for (int i = 0; i < l.size(); i++)
	{
		l[i].acceleration.y += g;
	}
}
void calculateGravityAttraction(Ball &b, std::vector<Ball> &l)
{
	for (int i = 0; i < l.size(); i++)
	{
		if (b.body.getPosition() != l[i].body.getPosition())
		{

			double dx = b.body.getPosition().x - l[i].body.getPosition().x;
			double dy = b.body.getPosition().y - l[i].body.getPosition().y;

			double r2 = dx*dx + dy*dy;
			double r = sqrt(r2);

			sf::Vector2f normDist;
			normDist.x = dx / r;
			normDist.y = dy / r;

			double force = G * l[i].mass() / r2;

			b.acceleration.x -= force * normDist.x;
			b.acceleration.y -= force * normDist.y;
		}
	}
}
void calculateGravityWellAttraction(std::vector<Ball> &l, std::vector<gravityWell> &gl)
{
	for (int q = 0; q < gl.size(); q++)
	{
		for (int i = 0; i < l.size(); i++)
		{
			double dx = l[i].body.getPosition().x - gl[q].pos.x;
			double dy = l[i].body.getPosition().y - gl[q].pos.y;
			double r2 = dx*dx + dy*dy;
			double r = sqrt(r2);

			if (r > GRAVITY_WELL_CUTOFF_CLOSE_DIST && r < gl[q].rad)
			{
				double acc = G * gl[q].strength / r;

				sf::Vector2f normDist;
				normDist.x = dx / r;
				normDist.y = dy / r;

				l[i].acceleration.x -= acc * normDist.x;
				l[i].acceleration.y -= acc * normDist.y;
			}
		}
	}
}
void calcForces(std::vector<Ball> &l)
{
	for (int i = 0; i < l.size(); i++) l[i].acceleration = sf::Vector2f(0, 0);
	if (GRAVITY_VERT_ENABLE) gravityVerticalBalls(ballList);
	if (GRAVITY_BETWEEN_BALLS)
	{
		for (int i = 0; i < l.size(); i++)
		{
			calculateGravityAttraction(l[i], l);
		}
	}
	calculateGravityWellAttraction(ballList, gravityWellList);
}
void updateVelocities(std::vector<Ball> &l)
{
	for (int i = 0; i < l.size(); i++)
	{
		l[i].velocity.x += l[i].acceleration.x*dt;
		l[i].velocity.y += l[i].acceleration.y*dt;

		l[i].acceleration.x = 0;
		l[i].acceleration.y = 0;
	}
}
void moveBalls(std::vector<Ball> &l)
{
	for (int i = 0; i < l.size(); i++)
	{
		sf::Vector2f dPos;
		dPos.x = l[i].velocity.x * dt;
		dPos.y = l[i].velocity.y * dt;

		l[i].body.move(dPos);
	}
}
void drawBalls(std::vector<Ball> &l, sf::RenderWindow &w)
{
	for (int i = 0; i < l.size(); i++)
	{
		w.draw(l[i].body);
	}
}
void addBallWithMouse(std::vector<Ball> &l, sf::RenderWindow &w, sf::Vector2f vel, sf::Vector2f pos)
{
	Ball newBall(pos, newBallRad);
	newBall.velocity = vel;
	l.push_back(newBall);
}
void addBallDrawMouseVector(sf::RenderWindow &w)
{
	if (newObjectInfoVector.z == 0 || !sf::Mouse::isButtonPressed(sf::Mouse::Left) || (!sf::Keyboard::isKeyPressed(sf::Keyboard::Q) && !sf::Keyboard::isKeyPressed(sf::Keyboard::E))) return;

	sf::Vector2f mPos = w.mapPixelToCoords(sf::Mouse::getPosition(w), w.getView());

	//CHOOSING COLOR BASED ON FUNCTION
	sf::Color col;
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::E)) col = sf::Color(0, 191, 255);
	else col = sf::Color::Red;

	//SLINGSHOT VECTOR
	sf::Vertex vector[2]
	{
		sf::Vertex(sf::Vector2f(newObjectInfoVector.x,newObjectInfoVector.y),col),
		sf::Vertex(mPos,col)
	};
	w.draw(vector, 2, sf::Lines);

	//BALL
	sf::CircleShape predBall(newBallRad);
	predBall.setRadius(newBallRad);
	predBall.setOrigin(newBallRad, newBallRad);
	predBall.setOutlineThickness(-newBallRad / 6.0);
	predBall.setFillColor(sf::Color(0, 0, 0, 0));
	predBall.setOutlineColor(col);
	predBall.setPointCount(200);
	predBall.setPosition(sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y));
	w.draw(predBall);

	//PREDICTION
	sf::Vector2f predictionPos;
	predictionPos.x = newObjectInfoVector.x;
	predictionPos.y = newObjectInfoVector.y;
	sf::Vector2f predictionVel;
	predictionVel.x = -slingShotStrength*(mPos.x - newObjectInfoVector.x);
	predictionVel.y = -slingShotStrength*(mPos.y - newObjectInfoVector.y);
	double predAlpha = PRED_VECTORS_START_ALPHA;
	for (int i = 0; i < PRED_VECTORS_AMOUNT; i++)
	{

		sf::Vertex vector[2]
		{
			sf::Vertex(predictionPos,sf::Color(col.r,col.g,col.b,predAlpha)),
			sf::Vertex(sf::Vector2f(predictionPos.x + predictionVel.x*PRED_VECTOR_TIME_STEP,predictionPos.y + predictionVel.y*PRED_VECTOR_TIME_STEP),sf::Color(col.r,col.g,col.b,predAlpha - PRED_VECTORS_START_ALPHA / PRED_VECTORS_AMOUNT))
		};
		w.draw(vector, 2, sf::Lines);

		predictionPos.x += predictionVel.x*PRED_VECTOR_TIME_STEP;
		predictionPos.y += predictionVel.y*PRED_VECTOR_TIME_STEP;
		if (GRAVITY_VERT_ENABLE) predictionVel.y += g*PRED_VECTOR_TIME_STEP;
		predAlpha -= PRED_VECTORS_START_ALPHA / PRED_VECTORS_AMOUNT;



	}
}
void checkIfBallCollidesFurther(Ball & b);

/*BALL COLLISION*/
bool isCollisionBetweenBalls(Ball &b1, Ball &b2)
{
	double dx = b1.body.getPosition().x - b2.body.getPosition().x;
	double dy = b1.body.getPosition().y - b2.body.getPosition().y;
	double dist2 = dx*dx + dy*dy;

	double sumRadi = b1.body.getRadius() + b2.body.getRadius();

	if (sumRadi*sumRadi < (dist2-1)) return false;
	else return true;
}
void collideBalls(Ball &b1, Ball &b2)
{
	//VECTOR BETWEEN BALLS
	sf::Vector2f normalVector = b1.body.getPosition() - b2.body.getPosition();
	double dist = sqrt(normalVector.x*normalVector.x + normalVector.y*normalVector.y);
	if (dist < 0.1) dist = 0.1;
	sf::Vector2f normalVectorNormalized = sf::Vector2f(normalVector.x / dist, normalVector.y /dist);

	//VELOCITY VECTOR
	sf::Vector2f deltaVelVector = b1.velocity - b2.velocity;
	b1.velocity.x = ball_bouncyness*(b1.velocity.x * (b1.mass() - b2.mass()) + (2 * b2.mass() * b2.velocity.x)) / (b1.mass() + b2.mass());
	b1.velocity.y = ball_bouncyness*(b1.velocity.y * (b1.mass() - b2.mass()) + (2 * b2.mass() * b2.velocity.y)) / (b1.mass() + b2.mass());
	b2.velocity.x = -ball_bouncyness*(b2.velocity.x * (b2.mass() - b1.mass()) + (2 * b1.mass() * b1.velocity.x)) / (b1.mass() + b2.mass());
	b2.velocity.y = -ball_bouncyness*(b2.velocity.y * (b2.mass() - b1.mass()) + (2 * b1.mass() * b1.velocity.y)) / (b1.mass() + b2.mass());

	//MOVING BALLS AWAY
	double travelDist = (collision_extra_dist + (b1.body.getRadius() + b2.body.getRadius())) - dist;
	b1.body.move(sf::Vector2f(0.5*travelDist * normalVectorNormalized.x, 0.5*travelDist * normalVectorNormalized.y));
	b2.body.move(sf::Vector2f(-0.5*travelDist * normalVectorNormalized.x, -0.5*travelDist * normalVectorNormalized.y));
}
void checkBallToBallCollision()
{

	for (int i = 0; i < ballList.size(); i++)
	{
		for (int j = 0; j < ballList.size(); j++)
		{
			if (isCollisionBetweenBalls(ballList[i], ballList[j]) && i != j)
			{
				collideBalls(ballList[i], ballList[j]);
				checkIfBallCollidesFurther(ballList[i]);
				checkIfBallCollidesFurther(ballList[j]);
				calcForces(ballList);
			}
		}
	}
}
void checkIfBallCollidesFurther(Ball & b)
{
	for (int i = 0; i < ballList.size(); i++)
	{
		if (isCollisionBetweenBalls(b, ballList[i]) && b.body.getPosition() != ballList[i].body.getPosition())
		{
			collideBalls(b, ballList[i]);
			checkIfBallCollidesFurther(b);
			checkIfBallCollidesFurther(ballList[i]);
		}
	}
}

/*WALLS*/
void drawWalls(std::vector<Wall> &l, sf::RenderWindow &w)
{
	for (int i = 0; i < l.size(); i++)
	{
		sf::Color col;
		if (l[i].killer) col = sf::Color(255, 0, 0,170);
		else col = sf::Color::White;


		sf::Vertex line[2]
		{
			sf::Vertex(l[i].pos,col),
			sf::Vertex(sf::Vector2f(l[i].pos.x + cos(l[i].angle) * l[i].length/2,l[i].pos.y + sin(l[i].angle) * l[i].length/2),col)
		};
		w.draw(line, 2, sf::Lines);

		sf::Vertex line2[2]
		{
			sf::Vertex(l[i].pos,col),
			sf::Vertex(sf::Vector2f(l[i].pos.x - cos(l[i].angle) * l[i].length/2,l[i].pos.y - sin(l[i].angle) * l[i].length/2),col)
		};
		w.draw(line2, 2, sf::Lines);
	}
}
bool isCollisionBetweenBallAndWall(Ball &b, Wall &w)
{

	//PROJECTING DELTAPOSITIONVECTOR ONTO WALL
	sf::Vector2f y = b.body.getPosition() - w.pos;
	sf::Vector2f wallNorm = sf::Vector2f(cos(w.angle+PI/2), sin(w.angle+PI/2));
	sf::Vector2f u = sf::Vector2f(cos(w.angle), sin(w.angle));
	double u_dot_y = y.x*u.x + y.y*u.y;
	sf::Vector2f y_proj_u = sf::Vector2f(u_dot_y * u.x, u_dot_y * u.y);
	
	if (sqrt(y_proj_u.x*y_proj_u.x + y_proj_u.y*y_proj_u.y) > w.length/2+b.body.getRadius()) return false;

	sf::Vector2f y_hat = y - y_proj_u;

	if (sqrt(y_hat.x*y_hat.x + y_hat.y*y_hat.y) > b.body.getRadius()) return false;

	//MOVING BALL SO IT CAN'T GET STUCK
	if (wallNorm.x*y.x+ wallNorm.y*y.y <= 0) b.body.setPosition(w.pos + y_proj_u - sf::Vector2f(cos(w.angle + PI*0.5)*b.body.getRadius(), sin(w.angle + PI*0.5)*b.body.getRadius()));
	else b.body.setPosition(w.pos + y_proj_u + sf::Vector2f(cos(w.angle + PI*0.5)*b.body.getRadius(), sin(w.angle + PI*0.5)*b.body.getRadius()));


	return true;


}
void collideBallWithWall(Ball &b, Wall &w)
{
	sf::Vector2f ballVel = b.velocity;
	sf::Vector2f bVel_proj_wall_norm;
	sf::Vector2f bVel_proj_wall_tang;

	//PROJECTING VELOCITYVECTOR ONTO WALL NORMAL VECTOR

	sf::Vector2f u = sf::Vector2f(cos(w.angle + PI / 2), sin(w.angle + PI / 2));
	double u_dot_y = ballVel.x*u.x + ballVel.y*u.y;
	bVel_proj_wall_norm = sf::Vector2f(wall_bouncyness*u_dot_y * u.x, wall_bouncyness*u_dot_y * u.y);

	//PROJECTING VELOCITYVECTOR ONTO WALL TANGENTVECTOR
	u = sf::Vector2f(cos(w.angle), sin(w.angle));
	u_dot_y = ballVel.x*u.x + ballVel.y*u.y;
	bVel_proj_wall_tang = sf::Vector2f(u_dot_y * u.x, u_dot_y * u.y);

	//SETTING NEW VELOCITY
	b.velocity = bVel_proj_wall_tang - bVel_proj_wall_norm;

}
void checkBallToWallCollision()
{
	for (int i = 0; i < ballList.size(); i++)
	{
		for (int p = 0; p < wallList.size(); p++)
		{
			if (isCollisionBetweenBallAndWall(ballList[i], wallList[p]))
			{
				if (wallList[p].killer)
				{
					auto it = ballList.begin() + i;
					ballList.erase(it);
					break;
				}
				else
				{
					collideBallWithWall(ballList[i], wallList[p]);
				}
			}
		}
	}
}
void addWallWithMouse(std::vector<Wall> &l, sf::RenderWindow &w, sf::Vector2f pos, bool kill)
{
	sf::Vector2f mPos = w.mapPixelToCoords(sf::Mouse::getPosition(w), w.getView());
	sf::Vector2f deltaPos = mPos - pos;
	double angle = atan2(deltaPos.y, deltaPos.x);
	double len = 2*sqrt(deltaPos.x*deltaPos.x + deltaPos.y*deltaPos.y);

	l.push_back(Wall(pos, len, angle, kill));
}
void addWallDrawNew(sf::RenderWindow &w)
{

	if (sf::Keyboard::isKeyPressed(sf::Keyboard::T) && !sf::Mouse::isButtonPressed(sf::Mouse::Left) && lastMouseFreeHand != sf::Vector2f(0,0))
	{

		sf::Vector2f loggedPos = sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y);
		sf::Vector2f mPos = w.mapPixelToCoords(sf::Mouse::getPosition(w), w.getView());
		sf::Vector2f deltaPos = mPos - loggedPos;

		sf::Vertex vizualiser2[2]
		{
			sf::Vertex(loggedPos,sf::Color(255, 165, 0)),
			sf::Vertex(mPos,sf::Color(255, 165, 0))
		};
		w.draw(vizualiser2, 2, sf::Lines);
	}

	if (newObjectInfoVector.z == 0 || !sf::Mouse::isButtonPressed(sf::Mouse::Left) || (!sf::Keyboard::isKeyPressed(sf::Keyboard::W) && !sf::Keyboard::isKeyPressed(sf::Keyboard::R) && !sf::Keyboard::isKeyPressed(sf::Keyboard::T))) return;
	
	sf::Color col;
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::W) || sf::Keyboard::isKeyPressed(sf::Keyboard::T)) col = sf::Color(255, 165, 0);
	else col = sf::Color(150, 0, 0, 150);

	if (sf::Keyboard::isKeyPressed(sf::Keyboard::W) || sf::Keyboard::isKeyPressed(sf::Keyboard::R))
	{

		sf::Vector2f loggedPos = sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y);
		sf::Vector2f mPos = w.mapPixelToCoords(sf::Mouse::getPosition(w), w.getView());
		sf::Vector2f deltaPos = mPos - loggedPos;

		sf::Vertex vizualiser[2]
		{
			sf::Vertex(loggedPos - deltaPos,col),
			sf::Vertex(loggedPos,col)
		};
		w.draw(vizualiser, 2, sf::Lines);

		sf::Vertex vizualiser2[2]
		{
			sf::Vertex(loggedPos + deltaPos,col),
			sf::Vertex(loggedPos,col)
		};
		w.draw(vizualiser2, 2, sf::Lines);
	}
	
}

/*SPAWNERS*/
void drawSpawners(std::vector<ballSpawner> &l, sf::RenderWindow &w)
{
	for (int i = 0; i < l.size(); i++)
	{
		sf::Vertex viz[2]
		{
			sf::Vertex(l[i].pos,sf::Color(0,191,255)),
			sf::Vertex(l[i].pos + sf::Vector2f(l[i].spawnVel.x+1,l[i].spawnVel.y+1),sf::Color(0,191,255,0))
		};
		w.draw(viz, 2, sf::Lines);
	}
}
void updateSpawners(std::vector<ballSpawner> &l, std::vector<Ball> &bl)
{
	for (int i = 0; i < l.size(); i++)
	{
		l[i].timeSinceLastSpawn += dt;
		if (l[i].timeSinceLastSpawn > l[i].delay)
		{
			l[i].timeSinceLastSpawn = 0;
			Ball newBall(l[i].pos, l[i].sizeBallSpawn);
			newBall.velocity = l[i].spawnVel;

			bl.push_back(newBall);
		}
	}
}
void addSpawnerWithMouse(std::vector<ballSpawner> &l, sf::RenderWindow &w, sf::Vector2f pos)
{
	sf::Vector2f mPos = w.mapPixelToCoords(sf::Mouse::getPosition(w), w.getView());
	sf::Vector2f deltaPos = pos - mPos;


	l.push_back(ballSpawner(pos,sf::Vector2f(deltaPos.x*slingShotStrength,deltaPos.y*slingShotStrength), newBallRad, 10));
}

/*GRAVITYWELL*/
void drawGravityWells(std::vector<gravityWell> &l, sf::RenderWindow &w)
{

	//WELLS
	for (int i = 0; i < l.size(); i++)
	{

		double angle;
		double dist = GRAVITY_WELL_CUTOFF_CLOSE_DIST;
		double acc = G*l[i].strength / dist;
		int numberLines = GRAVITY_FIELD_LINE_START_NUMBER;

		//LAYERS
		while (dist < l[i].rad)
		{
			angle = 0;

			//LINES
			for (int e = 0; e < numberLines; e++)
			{
				sf::Vertex viz[2]
				{
					sf::Vertex(l[i].pos + sf::Vector2f(dist*cos(angle),dist*sin(angle)),sf::Color(50,255,148,170)),
					sf::Vertex(l[i].pos + sf::Vector2f((dist - 0.1*acc)*cos(angle),(dist - 0.1*acc)*sin(angle)),sf::Color(50,255,148,0))
				};
				w.draw(viz, 2, sf::Lines);

				angle += 2 * PI / numberLines;
			}

			dist += 15;
			numberLines += 3;
			acc = g*l[i].strength / dist;

		}

		//CIRCLE
		sf::CircleShape outerCircle(l[i].rad);
		outerCircle.setOrigin(l[i].rad, l[i].rad);
		outerCircle.setPosition(l[i].pos);
		outerCircle.setFillColor(sf::Color(0, 0, 0, 0));
		outerCircle.setOutlineThickness(-1);
		outerCircle.setOutlineColor(sf::Color(50, 255, 148, 70));
		w.draw(outerCircle);
	}
}
void addGravityWellWithMouse(std::vector<gravityWell> &l, sf::RenderWindow &w, sf::Vector2f pos)
{
	sf::Vector2f mPos = w.mapPixelToCoords(sf::Mouse::getPosition(w), w.getView());
	sf::Vector2f deltaPos = pos - mPos;
	double r = sqrt(deltaPos.x*deltaPos.x + deltaPos.y*deltaPos.y);

	l.push_back(gravityWell(pos, 200 * newBallRad, r));
}
void addGravityWellDrawNew(sf::RenderWindow &w)
{
	if (newObjectInfoVector.z == 0 ) return;


	if (sf::Keyboard::isKeyPressed(sf::Keyboard::Y))
	{

		sf::Vector2f loggedPos = sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y);
		sf::Vector2f mPos = w.mapPixelToCoords(sf::Mouse::getPosition(w), w.getView());
		sf::Vector2f deltaPos = mPos - loggedPos;
		double r = sqrt(deltaPos.x*deltaPos.x + deltaPos.y*deltaPos.y);

		sf::CircleShape outerCircle(r);
		outerCircle.setOrigin(r,r);
		outerCircle.setPosition(loggedPos);
		outerCircle.setFillColor(sf::Color(0, 0, 0, 0));
		outerCircle.setOutlineThickness(-1);
		outerCircle.setOutlineColor(sf::Color(50, 255, 148, 100));
		w.draw(outerCircle);
	}

}


/*ICONS*/
void drawIcons(sf::RenderWindow &w)
{

	if (GRAVITY_VERT_ENABLE)
	{
		sf::Vertex arrow1[2]
		{
			sf::Vertex(sf::Vector2f(30,26),sf::Color::Red),
			sf::Vertex(sf::Vector2f(30,36),sf::Color::Red)
		};

		sf::Vertex arrow2[2]
		{
			sf::Vertex(sf::Vector2f(33,33),sf::Color::Red),
			sf::Vertex(sf::Vector2f(30,36),sf::Color::Red)
		};

		sf::Vertex arrow3[2]
		{
			sf::Vertex(sf::Vector2f(27,33),sf::Color::Red),
			sf::Vertex(sf::Vector2f(30,36),sf::Color::Red)
		};
		
		w.draw(arrow1, 2,sf::Lines);
		w.draw(arrow2, 2, sf::Lines);
		w.draw(arrow3, 2, sf::Lines);

	}

	if (GRAVITY_BETWEEN_BALLS)
	{

		sf::Vertex arrow[2]
		{
			sf::Vertex(sf::Vector2f(42,31),sf::Color::Red),
			sf::Vertex(sf::Vector2f(52,31),sf::Color::Red)
		};

		w.draw(arrow, 2, sf::Lines);

		sf::CircleShape circ1(3);
		circ1.setFillColor(sf::Color::Yellow);

		circ1.setPosition(40, 29);
		w.draw(circ1);

		circ1.setPosition(50, 29);
		w.draw(circ1);
	}

	if (PAUSE)
	{
		sf::Vertex arrow[2]
		{
			sf::Vertex(sf::Vector2f(64,28),sf::Color::Green),
			sf::Vertex(sf::Vector2f(64,35),sf::Color::Green)
		};

		sf::Vertex arrow2[2]
		{
			sf::Vertex(sf::Vector2f(68,28),sf::Color::Green),
			sf::Vertex(sf::Vector2f(68,35),sf::Color::Green)
		};

		w.draw(arrow, 2, sf::Lines);
		w.draw(arrow2, 2, sf::Lines);
	}
}

//WINDOW
sf::RenderWindow window(sf::VideoMode(1920, 1080), "Balls", sf::Style::Fullscreen);
sf::Event event;

//START
int main()
{
	window.setFramerateLimit(60);

	//SIDEWALLS
	wallList.push_back(Wall(sf::Vector2f(window.getSize().x/2, window.getSize().y-5), window.getSize().x, 0,false));
	wallList.push_back(Wall(sf::Vector2f(window.getSize().x / 2, 5), window.getSize().x, 0,false));
	wallList.push_back(Wall(sf::Vector2f(5, window.getSize().y/2), window.getSize().y, 0.5*PI,false));
	wallList.push_back(Wall(sf::Vector2f(window.getSize().x - 5, window.getSize().y / 2), window.getSize().y, 0.5*PI,false));

	while (window.isOpen())
	{

		while (window.pollEvent(event))
		{

			//EXIT GAME
			if (event.type == sf::Event::Closed) window.close();

			//ADD OBJECT WITH MOUSE
			if (event.type == sf::Event::MouseButtonPressed)
			{
				if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
				{

					newObjectInfoVector.x = window.mapPixelToCoords(sf::Mouse::getPosition(window), window.getView()).x;
					newObjectInfoVector.y = window.mapPixelToCoords(sf::Mouse::getPosition(window), window.getView()).y;
					newObjectInfoVector.z = 1;

					if (sf::Keyboard::isKeyPressed(sf::Keyboard::T) && lastMouseFreeHand != sf::Vector2f(0,0))
					{
						sf::Vector2f mPos = window.mapPixelToCoords(sf::Mouse::getPosition(window), window.getView());
						sf::Vector2f dPos = mPos - lastMouseFreeHand;

						Wall newWall(lastMouseFreeHand+sf::Vector2f(dPos.x/2,dPos.y/2),sqrt(dPos.x*dPos.x+dPos.y*dPos.y),atan2(dPos.y,dPos.x),false);
						wallList.push_back(newWall);
						
						lastMouseFreeHand = mPos;
					}
					else if (sf::Keyboard::isKeyPressed(sf::Keyboard::T))
					{
						sf::Vector2f mPos = window.mapPixelToCoords(sf::Mouse::getPosition(window), window.getView());
						lastMouseFreeHand = mPos;
					}
				}
			}
			if (event.type == sf::Event::MouseButtonReleased)
			{
				if (!sf::Mouse::isButtonPressed(sf::Mouse::Left) && newObjectInfoVector.z != 0 && sf::Keyboard::isKeyPressed(sf::Keyboard::Q))
				{
					sf::Vector2f mouseDelta;
					mouseDelta.x = -slingShotStrength*(window.mapPixelToCoords(sf::Mouse::getPosition(window), window.getView()).x - newObjectInfoVector.x);
					mouseDelta.y = -slingShotStrength*(window.mapPixelToCoords(sf::Mouse::getPosition(window), window.getView()).y - newObjectInfoVector.y);
					addBallWithMouse(ballList, window, mouseDelta,sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y));
					newObjectInfoVector.z = 0;
				}

				if (!sf::Mouse::isButtonPressed(sf::Mouse::Left) && newObjectInfoVector.z != 0 && sf::Keyboard::isKeyPressed(sf::Keyboard::W))
				{
					addWallWithMouse(wallList, window, sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y), false);
					newObjectInfoVector.z = 0;
				}

				if (!sf::Mouse::isButtonPressed(sf::Mouse::Left) && newObjectInfoVector.z != 0 && sf::Keyboard::isKeyPressed(sf::Keyboard::E))
				{
					addSpawnerWithMouse(spawnerList, window, sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y));
					newObjectInfoVector.z = 0;
				}

				if (!sf::Mouse::isButtonPressed(sf::Mouse::Left) && newObjectInfoVector.z != 0 && sf::Keyboard::isKeyPressed(sf::Keyboard::R))
				{
					addWallWithMouse(wallList, window, sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y), true);
					newObjectInfoVector.z = 0;
				}

				if (!sf::Mouse::isButtonPressed(sf::Mouse::Left) && newObjectInfoVector.z != 0 && sf::Keyboard::isKeyPressed(sf::Keyboard::Y))
				{
					addGravityWellWithMouse(gravityWellList, window, sf::Vector2f(newObjectInfoVector.x, newObjectInfoVector.y));
					newObjectInfoVector.z = 0;
				}

			}

			//SCROLLING CHANG0ES SIZE OF BALL && STRENGTH OF GRAVITYFIELD
			if (event.type == sf::Event::MouseWheelMoved)
			{
				double mouseDelta = event.mouseWheel.delta;

				if (mouseDelta > 0 && newBallRad < BALL_RAD_MAX) newBallRad += mouseDelta;
				else if (mouseDelta < 0 && newBallRad > BALL_RAD_MIN) newBallRad += mouseDelta;

				if (newBallRad > BALL_RAD_MAX) newBallRad = BALL_RAD_MAX;
				if (newBallRad < BALL_RAD_MIN) newBallRad = BALL_RAD_MIN;


			}

			//PAUSE
			if (event.type == sf::Event::KeyPressed && sf::Keyboard::isKeyPressed(sf::Keyboard::P))
			{
				if (!PAUSE)PAUSE = true;
				else PAUSE = false;
			}

			//CLEAR
			if (event.type == sf::Event::KeyPressed && sf::Keyboard::isKeyPressed(sf::Keyboard::C))
			{
				if (!sf::Keyboard::isKeyPressed(sf::Keyboard::Q)) ballList.clear();
				if (!sf::Keyboard::isKeyPressed(sf::Keyboard::W)) wallList.clear();
				if (!sf::Keyboard::isKeyPressed(sf::Keyboard::E)) spawnerList.clear();
				if (!sf::Keyboard::isKeyPressed(sf::Keyboard::Y)) gravityWellList.clear();


				//SIDEWALLS
				wallList.push_back(Wall(sf::Vector2f(window.getSize().x / 2, window.getSize().y - 5), window.getSize().x, 0, false));
				wallList.push_back(Wall(sf::Vector2f(window.getSize().x / 2, 5), window.getSize().x, 0, false));
				wallList.push_back(Wall(sf::Vector2f(5, window.getSize().y / 2), window.getSize().y, 0.5*PI, false));
				wallList.push_back(Wall(sf::Vector2f(window.getSize().x - 5, window.getSize().y / 2), window.getSize().y, 0.5*PI, false));

			}

			//GRAVITY SWITCHES
			if (event.type == sf::Event::KeyPressed && sf::Keyboard::isKeyPressed(sf::Keyboard::B))
			{
				if (GRAVITY_BETWEEN_BALLS) GRAVITY_BETWEEN_BALLS = false;
				else GRAVITY_BETWEEN_BALLS = true;
			}
			if (event.type == sf::Event::KeyPressed && sf::Keyboard::isKeyPressed(sf::Keyboard::V))
			{
				if (GRAVITY_VERT_ENABLE) GRAVITY_VERT_ENABLE = false;
				else GRAVITY_VERT_ENABLE = true;
			}
		}

		//GAME
		if (!PAUSE)
		{
			updateSpawners(spawnerList, ballList);
			checkBallToBallCollision();
			checkBallToWallCollision();
			calcForces(ballList);
			updateVelocities(ballList);
			moveBalls(ballList);
		}

		//CONTROLS
		if (!sf::Keyboard::isKeyPressed(sf::Keyboard::T)) lastMouseFreeHand = sf::Vector2f(0, 0);

		//DRAW TO WINDOW
		window.clear();
		drawBalls(ballList, window);
		drawWalls(wallList, window);
		drawSpawners(spawnerList, window);
		drawGravityWells(gravityWellList, window);
		drawIcons(window);
		addBallDrawMouseVector(window);
		addWallDrawNew(window);
		addGravityWellDrawNew(window);
		window.display();
	}

}