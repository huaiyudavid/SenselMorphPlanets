SenselDevice sensel;
boolean sensel_sensor_opened = false;
ArrayList<Planet> system = new ArrayList<Planet>();
int WINDOW_WIDTH = 1400;
int WINDOW_HEIGHT;

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
  noStroke();
}

float[] force = new float[1000];
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
      force[i] = (force[i] * 0.95f) + (contacts[i].total_force * 0.05f);
      float r = force[i] / 100f;
      int pX = (int) ((contacts[i].x_pos_mm / sensel.getSensorWidthMM())  * WINDOW_WIDTH);
      int pY = (int) ((contacts[i].y_pos_mm / sensel.getSensorHeightMM()) * WINDOW_HEIGHT);
      if (contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_END){
        system.add(new Planet(pX, pY, 0, 0, 0, 0, r));
      } else if ( contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_MOVE || contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_START){
        ellipse((float)pX, (float)pY, r, r);
      }
    }
  }
  for (Planet p : system){
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
