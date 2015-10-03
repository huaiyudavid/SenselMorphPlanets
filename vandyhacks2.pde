SenselDevice sensel;
boolean sensel_sensor_opened = false;
ArrayList<Planet> system = new ArrayList<Planet>();
int WINDOW_WIDTH = 1400;
int WINDOW_HEIGHT;
int FPS = 60;
  
void setup()
{
  //sensel setup
  DisposeHandler dh = new DisposeHandler(this);
  sensel = new SenselDevice(this);
  sensel_sensor_opened = sensel.openConnection();
  if(!sensel_sensor_opened)
  {
    println("Unable to open Sensel sensor!");
    exit();
    return;
  }
  sensel.setFrameContentControl(SenselDevice.SENSEL_FRAME_CONTACTS_FLAG);
  sensel.startScanning();
  //end sensel setup
  WINDOW_HEIGHT = (int) ((sensel.getSensorHeightMM() / sensel.getSensorWidthMM()) * 1400);
  size(WINDOW_WIDTH,WINDOW_HEIGHT, P3D);
  frameRate(FPS);
  noStroke();
}

float[] force = new float[1000];
int[] initialX = new int[1000];
int[] initialY = new int[1000];

void draw()
{
  if(!sensel_sensor_opened) {
    return;
  }
  background(0);
  SenselContact[] contacts = sensel.readContacts();
  if(contacts != null)
  {
    for (int i = 0; i < contacts.length; i++){
      int id = contacts[i].id;
      force[id] = (force[id] * 0.95f) + (contacts[i].total_force * 0.05f);
      float r = force[id] / 100f;
      int pX = (int) ((contacts[i].x_pos_mm / sensel.getSensorWidthMM())  * WINDOW_WIDTH);
      int pY = (int) ((contacts[i].y_pos_mm / sensel.getSensorHeightMM()) * WINDOW_HEIGHT);
      if (contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_END){
        system.add(new Planet(initialX[id], initialY[id], 0, (pX - initialX[id]) / 4f, (pY - initialY[id]) / 5f, 0, r));
        force[id] = 0;
      } else if ( contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_MOVE || contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_START){
        if (r >= 10){
          Planet temp = new Planet(initialX[id],initialY[id], 0, 0,0,0,r);
          for (int a = 0; a < system.size(); a++) {
            if (temp.collides(system.get(a))){
              force[id] = (float)(distance(system.get(a).x, system.get(a).y, system.get(a).z, temp.x, temp.y, temp.z) - system.get(a).radius) * 100f;
            }
          }
          temp.radius = force[id] / 100f;
          temp.draw();
        }
        if (contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_START){
          initialX[id] = pX;
          initialY[id] = pY;
        }
      }
    }
  }
  
  for (int i = 0; i < system.size(); i++){
    Planet pl = system.get(i);
    if (pl.radius < 10 || ((pl.x + pl.radius < 0 || pl.x - pl.radius > WINDOW_WIDTH) && (pl.y + pl.radius < 0 || pl.y - pl.radius > WINDOW_HEIGHT))){
        system.remove(i);
        i--;
    }
  }
  
  for (int i = 0; i < system.size(); i++) {
    Planet p1 = system.get(i);
    for (int j = i+1; j < system.size(); j++) {
      Planet p2 = system.get(j);
      if (p1.collides(p2)) {
        elasticCollision(p1, p2);
      }
    }
  }
  
  for (Planet p : system){
    p.update(1.0 / FPS, system);
    p.draw();
  }
    
}

void keyPressed(){
  if (key=='q'||key=='Q'){
    exit();
    return;
  }
}

public class DisposeHandler 
{   
  DisposeHandler(PApplet pa)
  {
    pa.registerMethod("dispose", this);
  }  
  public void dispose()
  {      
    println("Closing sketch");
    if(sensel_sensor_opened)
    {
      sensel.stopScanning();
      sensel.closeConnection();
    }
  }
}
