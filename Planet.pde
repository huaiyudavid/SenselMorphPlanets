
public class Planet
{
  double x, y, z;
  double velX, velY, velZ;
  float radius, mass; //mass = radius
  PShape shape;
  
  public Planet(double x, double y, double z, double velX, double velY, double velZ, float r)
  {
    this.x = x;
    this.y = y;
    this.z = z;
    this.velX = velX;
    this.velY = velY;
    this.velZ = velZ;
    this.radius = r;
    this.mass = r;
    
    shape = createShape(ELLIPSE, 0, 0, radius, radius);
  }
  
  public void draw()
  {
    shape(shape, (float)x, (float)y);
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
        g.multiply((p.radius * this.radius) / (distance(p.x, p.y, this.x, this.y) * distance(p.x, p.y, this.x, this.y)));
        netForce.add(g);
      }
    }
  }
  
  public boolean collides(Planet other)
  {
    double distance = Planet.distance(this.x, this.y, this.z, other.x, other.y, other.z);
    double collisionDistance = (double)(radius + other.radius);
    if (distance <= collisionDistance)
      return true;
    else
      return false;
  }
  
  public static double distance(double x1, double y1, double z1, double x2, double y2, double z2)
  {
     return Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2) + Math.pow(z2-z1, 2));
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
    double magnitude = Planet.distance(i, j, k, 0, 0, 0);
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
