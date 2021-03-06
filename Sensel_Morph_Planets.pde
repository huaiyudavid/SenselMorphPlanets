SenselDevice sensel;
boolean sensel_sensor_opened = false;
ArrayList<Planet> system = new ArrayList<Planet>();
int WINDOW_WIDTH = 1400;
int WINDOW_HEIGHT = 730;
int FPS = 60;

void setup()
{
  size(WINDOW_WIDTH, WINDOW_HEIGHT, P3D);

  //drawing setup
  //ortho();
  //perspective(PI/2.0, width/height, ((height/2.0) / tan(PI*60.0/360.0))/10.0, ((height/2.0) / tan(PI*60.0/360.0))*10.0);  
  perspective();
  //frustum(0, width, height, 0, 10, 200);
    //sensel setup
    DisposeHandler dh = new DisposeHandler(this);
  sensel = new SenselDevice(this);
  sensel_sensor_opened = sensel.openConnection();
  if (!sensel_sensor_opened)
  {
    println("Unable to open Sensel sensor!");
    exit();
    return;
  }
  sensel.setFrameContentControl(SenselDevice.SENSEL_FRAME_CONTACTS_FLAG);
  sensel.startScanning();
  //end sensel setup
  frameRate(FPS);
  noStroke();
  
  float[] acc_data = sensel.readAccelerometerData();
  rotXi = acc_data[1];
  rotYi = acc_data[0];
  rotZi = 0;
  cameraPointi = new double[]{width/2, height/2, (height/2.0) / tan(PI*30.0 / 180.0)};
}

float[] force = new float[1000];
int[] initialX = new int[1000];
int[] initialY = new int[1000];
double rotXi, rotYi, rotZi;
double[] cameraPointi;
float[] acc_data_cache = null;

void draw()
{
  if (!sensel_sensor_opened) {
    return;
  }
  background(0);

  //setup drawing
  //directionalLight(255, 255, 255, 0, 0, -1);
  lights();
  sphereDetail(30);
  SenselContact[] contacts = sensel.readContacts();
  float[] acc_data = sensel.readAccelerometerData();
  if (acc_data_cache == null){
    acc_data_cache = acc_data;
  } else {
    for (int i = 0; i < acc_data.length; i++){
      acc_data_cache[i] = (acc_data_cache[i] * 0.95) + (acc_data[i] * 0.05);
    }
  }
  
  double rotXangle = (acc_data_cache[1] - rotXi) * PI/2;
  double rotYangle = (acc_data_cache[0] - rotZi) * PI/2;
  double rotZangle = 0;
  
  //double[] cameraP = cameraPointi;
  double[] cameraP = rotX(cameraPointi, rotXangle);
  //cameraP = rotY(cameraP, rotYangle);
  //println("Acc Data: (" + acc_data[0] + ", " + acc_data[1] + ", " + acc_data[2] + ")");
  
  camera((float)cameraP[0], (float)cameraP[1], (float)cameraP[2], width/2.0, height/2.0, 0, 0, 1, 0);

  if (contacts != null && acc_data != null)
  {
    for (int i = 0; i < contacts.length; i++) {
      int id = contacts[i].id;
      force[id] = (force[id] * 0.95f) + (contacts[i].total_force * 0.05f);
      float r = force[id] / 100f;
      int pX = (int) ((contacts[i].x_pos_mm / sensel.getSensorWidthMM())  * WINDOW_WIDTH);
      int pY = (int) ((contacts[i].y_pos_mm / sensel.getSensorHeightMM()) * WINDOW_HEIGHT);
      if (contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_END) {
        system.add(new Planet(initialX[id], initialY[id], 0, (pX - initialX[id]) / 4f, (pY - initialY[id]) / 5f, 0, r));
        force[id] = 0;
      } else if ( contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_MOVE || contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_START) {
        if (r >= 10) {
          Planet temp = new Planet(initialX[id], initialY[id], 0, 0, 0, 0, r);
          for (int a = 0; a < system.size (); a++) {
            if (temp.collides(system.get(a))) {
              force[id] = (float)(distance(system.get(a).x, system.get(a).y, system.get(a).z, temp.x, temp.y, temp.z) - system.get(a).radius) * 100f;
              r = force[id] / 100f;
            }
          }
          temp = new Planet(initialX[id], initialY[id], 0, 0, 0, 0, r);
          temp.draw();
        }
        if (contacts[i].type == SenselDevice.SENSEL_EVENT_CONTACT_START) {
          initialX[id] = pX;
          initialY[id] = pY;
        }
      }
    }
  }

  for (int i = 0; i < system.size (); i++) {
    Planet pl = system.get(i);
    if (pl.radius < 10 || ((pl.x + pl.radius < -10000 || pl.x - pl.radius > WINDOW_WIDTH * 11) || (pl.y + pl.radius < -10000 || pl.y - pl.radius > WINDOW_HEIGHT * 11))) {
      system.remove(i);
      i--;
    }
  }

  for (int i = 0; i < system.size (); i++) {
    Planet p1 = system.get(i);
    for (int j = i+1; j < system.size (); j++) {
      Planet p2 = system.get(j);
      if (p1.collides(p2)) {
        elasticCollision(p1, p2);
      }
    }
  }

  for (Planet p : system) {
    p.update(3.0 / FPS, system);
    p.draw();
  }
}

void keyPressed() {
  if (key=='q'||key=='Q') {
    exit();
    return;
  }
  if (key=='r'||key=='R'){
    system = new ArrayList<Planet>();
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
    if (sensel_sensor_opened)
    {
      sensel.stopScanning();
      sensel.closeConnection();
    }
  }
}

