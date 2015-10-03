
public class Planet
{
  float x, y, z;
  float velX, velY, velZ;
  float radius, mass; //mass = radius
  PShape shape;
  
  public Planet(float x, float y, float z, float velX, float velY, float velZ, float r)
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
    shape(shape, x, y);
  }
  
  public void update(float dt, ArrayList<Planet> planets)
  {
    x += 
  }
}
