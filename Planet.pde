public static double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2) + Math.pow(z2-z1, 2));
}  

public class Planet
{
  double x, y, z;
  double velX, velY, velZ;
  float radius, mass; //mass = radius
  PShape shape;
  color c;
  final int G = 100;
  
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
