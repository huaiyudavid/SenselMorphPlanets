SenselDevice sensel;
boolean sensel_sensor_opened = false;
ArrayList<Planet> system = new ArrayList<Planet>();
int WINDOW_WIDTH = 1400;
int WINDOW_HEIGHT;
int FPS = 30;

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
  size(WINDOW_WIDTH,WINDOW_HEIGHT);
  framerate(FPS);
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
      force[i] = (force[i] > contacts[i].total_force) ? (force[i]) : (contacts[i].total_force);
      float r = force[i] / 100f;
      int pX = (int) ((contacts[i].x_pos_mm / sensel.getSensorWidthMM())  * WINDOW_WIDTH);
      int pY = (int) ((contacts[i].y_pos_mm / sensel.getSensorHeightMM()) * WINDOW_HEIGHT);
      if (contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_END){
        system.add(new Planet(initialX[i], initialY[i], 0, (pX - initialX[i]) / 5f, (pY - initialY[i]) / 5f, 0, r));
        force[i] = -1;
      } else if ( contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_MOVE || contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_START){
        ellipse(initialX[i], initialY[i], r, r);
        if (contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_START){
          initialX[i] = pX;
          initialY[i] = pY;
        }
      }
    }
  }
  for (Planet p : system){
    p.update(1.0 / FPS, system);
    ellipse((float)p.x, (float)p.y, (float)p.radius, (float)p.radius);
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
