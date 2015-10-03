public static double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2) + Math.pow(z2-z1, 2));
}

//Precondition: p1 collides with p2
public static void elasticCollision(Planet p1, Planet p2)
{
  double m1 = p1.mass, m2 = p2.mass;
  double r1 = p1.radius, r2 = p2.radius;
  double x1 = p1.x, y1 = p1.y, z1 = p1.z;
  double x2 = p2.x, y2 = p2.y, z2 = p2.z;
  double vx1 = p1.velX, vy1 = p1.velY, vz1 = p1.velZ;
  double vx2 = p2.velX, vy2 = p2.velY, vz2 = p2.velZ;

  double r12 = r1 + r2;
  double m21 = m2 / m1;
  double x21 = x2 - x1, y21 =y2 - y1, z21 = z2 - z1; 
  double vx21 = vx2 - vx1, vy21 = vy2 - vy1, vz21 = vz2 - vz1;
  double vx_cm = (m1*vx1+m2*vx2)/(m1+m2), vy_cm = (m1*vy1+m2*vy2)/(m1+m2), vz_cm = (m1*vz1+m2*vz2)/(m1+m2);

  //     **** calculate relative distance and relative speed ***
  double d = Math.sqrt(x21*x21 + y21*y21 + z21*z21);
  double v = Math.sqrt(vx21*vx21 + vy21*vy21 + vz21*vz21);

  //     **** return if relative speed = 0 ****
  if (v == 0)
    return;

  //     **** shift coordinate system so that ball 1 is at the origin ***
  x2 = x21;
  y2 = y21;
  z2 = z21;

  //     **** boost coordinate system so that ball 2 is resting ***
  vx1 = -vx21;
  vy1 = -vy21;
  vz1 = -vz21;

  //     **** find the polar coordinates of the location of ball 2 ***
  double theta2 = Math.acos(z2/d);
  double phi2;
  if (x2 == 0 && y2 == 0) 
  {
    phi2 = 0;
  } else 
  {
    phi2 = Math.atan2(y2, x2);
  }
  double st = Math.sin(theta2);
  double ct = Math.cos(theta2);
  double sp = Math.sin(phi2);
  double cp = Math.cos(phi2);

  //     **** express the velocity vector of ball 1 in a rotated coordinate
  //          system where ball 2 lies on the z-axis ******
  double vx1r = ct*cp*vx1 + ct*sp*vy1 - st*vz1;
  double vy1r = cp*vy1 - sp*vx1;
  double vz1r = st*cp*vx1 + st*sp*vy1 + ct*vz1;
  double fvz1r = vz1r/v;
  if (fvz1r > 1) {
    fvz1r = 1;
  }   // fix for possible rounding errors
  else if (fvz1r < -1) {
    fvz1r = -1;
  } 
  double thetav = Math.acos(fvz1r);
  double phiv;
  if (vx1r == 0 && vy1r == 0) {
    phiv = 0;
  } else {
    phiv = Math.atan2(vy1r, vx1r);
  }

  //     **** calculate the normalized impact parameter ***
  double dr = d*sin(thetav) / r12;

  //     **** calculate impact angles when balls collide ***
  double alpha = Math.asin(-dr);
  double beta = phiv;
  double sbeta = Math.sin(beta);
  double cbeta = Math.cos(beta);

  //     **** reverse the coordinate shift ***
  x2 += x1;
  y2 += y1;
  z2 += z1;

  //  ***  update velocities ***

  double a = Math.tan(thetav+alpha);

  double dvz2 = 2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*(1+m21));

  vz2r = dvz2;
  vx2r = a*cbeta*dvz2;
  vy2r = a*sbeta*dvz2;
  vz1r = vz1r-m21*vz2r;
  vx1r = vx1r-m21*vx2r;
  vy1r = vy1r-m21*vy2r;

  //     **** rotate the velocity vectors back and add the initial velocity
  //           vector of ball 2 to retrieve the original coordinate system ****

  vx1 = ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2;
  vy1 = ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2;
  vz1 = ct*vz1r-st*vx1r               +vz2;
  vx2 = ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2;
  vy2 = ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2;
  vz2 = ct*vz2r-st*vx2r               +vz2;
  
  //     **** update the planet velocities with the calculated velocities ****
  p1.velX = vx1;
  p1.velY = vy1;
  p1.velZ = vz1;
  p2.velX = vx2;
  p2.velY = vy2;
  p2.velZ = vz2;
}

public class Planet
{
  double x, y, z;
  double velX, velY, velZ;
  double radius, mass; //mass = radius
  PShape shape;
  color c;
  final int G = 100;

  public Planet(double x, double y, double z, double velX, double velY, double velZ, double r)
  {
    this.x = x;
    this.y = y;
    this.z = z;
    this.velX = velX;
    this.velY = velY;
    this.velZ = velZ;
    this.radius = r;
    this.mass = r;
    int red = (int)(r < 150 ? ((r-100) / 50f) * 255 : 255);
    int green = (int)(r < 100 ? ((r-50) / 50f) * 255 : ( r < 150 ? 255 : (((200 - r)  / 50f) * 255)));
    int blue = (int)(r < 50 ? (r / 50f) * 255 : ( r < 100 ? 255 : (((150 - r)  / 50f) * 255)));
    c = color(red,green,blue);
  }

  public void draw()
  {
    fill(c);
    ellipse((float)x, (float)y, (float)radius, (float)radius);
  }

  public void update(double dt, ArrayList<Planet> planets)
  {
    x += velX * dt;
    y += velY * dt;
    z += velZ * dt;

    Vector netForce = new Vector(0, 0, 0);

    for (Planet p : planets)
    {
      if (this != p)
      {
        Vector g = new Vector(p.x - this.x, p.y - this.y, p.z - this.z);
        g.normalize();
        g.multiply((p.radius * this.radius) / (distance(p.x, p.y, 0, this.x, this.y, 0) * distance(p.x, p.y, 0, this.x, this.y, 0)));
        netForce.add(g);
      }
    }
    netForce.multiply(G);
    velX += netForce.i / this.radius;
    velY += netForce.j / this.radius;
    velZ += netForce.k / this.radius;
  }

  public boolean collides(Planet other)
  {
    double distance = distance(this.x, this.y, this.z, other.x, other.y, other.z);
    double collisionDistance = (double)(radius + other.radius);
    if (distance <= collisionDistance)
      return true;
    else
      return false;
  }
}

public class Vector
{
  double i, j, k;

  public Vector(double i, double j, double k)
  {
    this.i = i;
    this.j = j;
    this.k = k;
  }

  public void add(Vector other)
  {
    i += other.i;
    j += other.j;
    k += other.k;
  }

  public void normalize()
  {
    double magnitude = distance(i, j, k, 0, 0, 0);
    i /= magnitude;
    j /= magnitude;
    k /= magnitude;
  }

  public void multiply(double other)
  {
    i *= other;
    j *= other;
    k *= other;
  }

  public double dot(Vector other)
  {
    return i*other.i + j*other.j + k*other.k;
  }
}

