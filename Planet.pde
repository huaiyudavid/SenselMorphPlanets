
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
    
  }
  
  public boolean collides(Planet other)
  {
    return false;
  }
  
  private double distance(double x1, double y1, double x2, double y2)
  {
     return Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2));
  }
}
