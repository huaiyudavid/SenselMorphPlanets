
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
        
        
      }
    }
  }
  
  public boolean collides(Planet other)
  {
    double distance = distance(this.x, this.y, other.x, other.y);
    double collisionDistance = (double)(radius + other.radius);
    if (distance <= collisionDistance)
      return true;
    else
      return false;
  }
  
  private double distance(double x1, double y1, double x2, double y2)
  {
     return Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2));
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
    
  }
  
}
