import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import processing.serial.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class Sensel_Morph_Planets extends PApplet {

SenselDevice sensel;
boolean sensel_sensor_opened = false;
ArrayList<Planet> system = new ArrayList<Planet>();
int WINDOW_WIDTH = 1400;
int WINDOW_HEIGHT = 730;
int FPS = 60;

public void setup()
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
  cameraPointi = new double[]{width/2, height/2, (height/2.0f) / tan(PI*30.0f / 180.0f)};
}

float[] force = new float[1000];
int[] initialX = new int[1000];
int[] initialY = new int[1000];
double rotXi, rotYi, rotZi;
double[] cameraPointi;
float[] acc_data_cache = null;

public void draw()
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
      acc_data_cache[i] = (acc_data_cache[i] * 0.95f) + (acc_data[i] * 0.05f);
    }
  }
  
  double rotXangle = (acc_data_cache[1] - rotXi) * PI/2;
  double rotYangle = (acc_data_cache[0] - rotZi) * PI/2;
  double rotZangle = 0;
  
  //double[] cameraP = cameraPointi;
  double[] cameraP = rotX(cameraPointi, rotXangle);
  //cameraP = rotY(cameraP, rotYangle);
  //println("Acc Data: (" + acc_data[0] + ", " + acc_data[1] + ", " + acc_data[2] + ")");
  
  camera((float)cameraP[0], (float)cameraP[1], (float)cameraP[2], width/2.0f, height/2.0f, 0, 0, 1, 0);

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
    p.update(3.0f / FPS, system);
    p.draw();
  }
}

public void keyPressed() {
  if (key=='q'||key=='Q') {
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
    if (sensel_sensor_opened)
    {
      sensel.stopScanning();
      sensel.closeConnection();
    }
  }
}

public static double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2) + Math.pow(z2-z1, 2));
}

public static double[] rotX(double[] point, double angle)
{
  double x = point[0], y = point[1], z = point[2];
  double newY = y * Math.cos(angle) - z * Math.sin(angle);
  double newZ = y * Math.sin(angle) + z * Math.cos(angle);
  return new double[]{x, newY, newZ};
}

public static double[] rotY(double[] point, double radians)
{
  double x = point[0], y = point[1], z = point[2];
  double newZ = z * Math.cos(radians) - x * Math.sin(radians);
  double newX = z * Math.sin(radians) + x * Math.cos(radians);
  return new double[]{newX, y, newZ};
}

public static double[] rotZ(double[] point, double radians)
{
  double x = point[0], y = point[1], z = point[2];
  double newX = x * Math.cos(radians) - y * Math.sin(radians);
  double newY = x * Math.sin(radians) + y * Math.cos(radians);
  return new double[]{newX, newY, z};
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
  double dr = d*Math.sin(thetav) / r12;

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

  double vz2r = dvz2;
  double vx2r = a*cbeta*dvz2;
  double vy2r = a*sbeta*dvz2;
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
  //p1.velZ = vz1;
  p2.velX = vx2;
  p2.velY = vy2;
  //p2.velZ = vz2;
}

public class Planet
{
  double x, y, z;
  double velX, velY, velZ;
  double radius, mass; //mass = radius
  PShape shape;
  int col;
  final int G = 1000;

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
    System.out.println(red + ", " + green + ", " + blue);
    col = color(red, green, blue);
    //shape = createShape(ELLIPSE, 0f, 0f, (float)r * 2, (float)r * 2);
    shape = createShape(SPHERE, (float)r);
  }

  public void draw()
  {
    shape.translate((float)(x + radius), (float)(y+ radius), (float)(z + radius));
    noStroke();
    fill(col);
    shape(shape);
    shapeMode(CENTER);
    //shape.translate(-(float)(x + radius), -(float)(y+ radius), -(float)(z + radius));
    shape.resetMatrix();
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
    velX += netForce.i / this.mass;
    velY += netForce.j / this.mass;
    velZ += netForce.k / this.mass;
  }

  public boolean collides(Planet other)
  {
    double distance = distance(this.x, this.y, this.z, other.x, other.y, other.z);
    double collisionDistance = radius + other.radius;
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

/**************************************************************************
 * Copyright 2015 Sensel, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 **************************************************************************/

/* 
  Sensel API for Processing.

  -> Outside of all functions:
  -> SenselDevice sensel;

  -> Inside setup():
  -> sensel = new SenselDevice(this);
  -> sensel.openConnection(); //returns true if successful, false if error
  -> sensel.setFrameContentControl(SenselDevice.SENSEL_FRAME_CONTACTS_FLAG); //enables contact sending
  -> sensel.startScanning(); //Start scanning
  
  -> Inside draw():
  -> SenselContact[] contacts = sensel.readContacts();

  -> Inside DisposeHandler:
  -> sensel.stopScanning();
  -> sensel.closeConnection();
*/



public class SenselContact
{
  int total_force;
  int uid;
  float area_mm_sq; // area in square mm
  float x_pos_mm; // x position in mm                                                                                                                                         
  float y_pos_mm; // y position in mm                                                                                                                                           
  float dx_mm; // change in x from last frame                                                                                                                                                                     
  float dy_mm; // change in y from last frame                                                                                                                                                                     
  float orientation_degrees; // angle from -90 to 90 degrees                                                                                                                          
  float major_axis_mm; // length of the major axis                                                                                                                                             
  float minor_axis_mm; // length of the minor axis                                                                                                                                             
  int id;
  int type;
}

public class SenselDevice
{
  final static byte SENSEL_REG_MAGIC                = (byte)0x00;
  final static byte SENSEL_REG_MAGIC_LENGTH         = (byte)0x06;
  final static byte SENSEL_REG_SCAN_CONTENT_CONTROL = (byte)0x24;
  final static byte SENSEL_REG_SCAN_ENABLED         = (byte)0x25;
  final static byte SENSEL_REG_SCAN_READ_FRAME      = (byte)0x26;
  final static byte SENSEL_REG_CONTACTS_MAX_COUNT   = (byte)0x40;
  final static byte SENSEL_REG_ACCEL_X              = (byte)0x60;
  final static byte SENSEL_REG_ACCEL_Y              = (byte)0x62;
  final static byte SENSEL_REG_ACCEL_Z              = (byte)0x64;
  final static byte SENSEL_REG_LED_BRIGHTNESS       = (byte)0x80;

  final static byte SENSEL_FRAME_CONTACTS_FLAG      = (byte)0x04;

  final static byte SENSEL_PT_FRAME     = 1;
  final static byte SENSEL_PT_READ_ACK  = 6;
  final static byte SENSEL_PT_WRITE_ACK = 10;

  final static int SENSEL_EVENT_CONTACT_INVALID = 0;
  final static int SENSEL_EVENT_CONTACT_START   = 1;
  final static int SENSEL_EVENT_CONTACT_MOVE    = 2;
  final static int SENSEL_EVENT_CONTACT_END     = 3;

  final static byte SENSEL_BOARD_ADDR   = (byte)0x01;
  final static byte SENSEL_READ_HEADER  = (byte)(SENSEL_BOARD_ADDR | (1 << 7));
  final static byte SENSEL_WRITE_HEADER = SENSEL_BOARD_ADDR;

  private Serial serial_port;
  private int sensor_max_x;
  private int sensor_max_y;
  private float sensor_width_mm;
  private float sensor_height_mm;
  private int sensor_max_contacts;
  private float sensor_x_to_mm_factor;
  private float sensor_y_to_mm_factor;
  private float sensor_orientation_to_degrees_factor = 1.0f/256.0f;
  private float sensor_area_to_mm_sq_factor = 1.0f/4096.0f;
  private PApplet parent;

  public SenselDevice(PApplet p)
  {
    parent = p; 
  }

  private boolean _checkForMagic(Serial port)
  {
    port.write(SENSEL_READ_HEADER);
    port.write(SENSEL_REG_MAGIC);
    port.write(SENSEL_REG_MAGIC_LENGTH);
    
    delay(500);
    
    //1-byte packet type, 2-byte size of payload, Payload, 1-byte checksum
    int magic_response_size = 4 + SENSEL_REG_MAGIC_LENGTH;
    
    if(port.available() < magic_response_size)
    {
      println("Magic not found!");
    }
    else
    {
      //println("Bytes available: " + port.available());

      //Check ACK
      if(port.readChar() != SENSEL_PT_READ_ACK) //Packet ACK
      {
        println("READ ACK NOT FOUND IN MAGIC PACKET!");
        return false;
      }

      //Check 2-byte packet size
      int packet_size = _convertBytesTo16((byte)port.readChar(), (byte)port.readChar());
      if(packet_size != SENSEL_REG_MAGIC_LENGTH)
      {
        println("LENGTH MISMATCH IN MAGIC PACKET! (Expected " + SENSEL_REG_MAGIC_LENGTH + ", received " + packet_size + ")");
        return false;
      }
      
      String magic = "";
      int checksum_calculated = 0;
      for(int i = 0; i < SENSEL_REG_MAGIC_LENGTH; i++)
      {
        char c = port.readChar();
        magic += c;
        checksum_calculated += c; 
      }
      checksum_calculated &= (0xFF);
      
      //Verify checksum
      int checksum_received = (int)port.readChar();
      if(checksum_received != checksum_calculated)
      {
        println("CHECKSUM MISMATCH IN MAGIC PACKET! (calculated " + (int)checksum_calculated + ", received " + (int)checksum_received + ")");
        return false;
      }
      
      if(magic.equals("S3NS31"))
      {
        println("MAGIC FOUND!");
        return true;
      }
      else
      {
        println("Invalid magic: " + magic);
      }
    }  
    return false;
  }

  public boolean openConnection()
  {
    String[] serial_list = Serial.list();
    serial_port = null;
    
    for(int i = 0; i < serial_list.length; i++)
    {
      println("Opening " + serial_list[i]);
      Serial curr_port;
      
      try{
        curr_port = new Serial(parent, serial_list[i], 115200);
      }
      catch(Exception e)
      {
        continue;
      }
      
      //Flush port
      curr_port.clear();
      
      if(_checkForMagic(curr_port))
      {
        serial_port = curr_port;
        sensor_max_x =  (256 * (_readReg(0x10, 1)[0] - 1));
        sensor_max_y =  (256 * (_readReg(0x11, 1)[0] - 1));
        
        println("Sensor Max X = " + sensor_max_x);
        println("Sensor Max Y = " + sensor_max_y);
        
        int [] sensor_width_arr  = _readReg(0x14, 4);
        int [] sensor_height_arr = _readReg(0x18, 4);
        
        //Convert from um to mm
        sensor_width_mm  = ((float)_convertBytesTo32((byte)sensor_width_arr[0],  (byte)sensor_width_arr[1],  (byte)sensor_width_arr[2],  (byte)sensor_width_arr[3]))  / 1000.0f;
        sensor_height_mm = ((float)_convertBytesTo32((byte)sensor_height_arr[0], (byte)sensor_height_arr[1], (byte)sensor_height_arr[2], (byte)sensor_height_arr[3])) / 1000.0f;
        
        println("Sensor Width = "  + sensor_width_mm  + " mm");
        println("Sensor Height = " + sensor_height_mm + " mm");
        
        sensor_max_contacts = _readReg(SENSEL_REG_CONTACTS_MAX_COUNT, 1)[0];
        println("Sensor Max Contacts = " + sensor_max_contacts);
        
        sensor_x_to_mm_factor = sensor_width_mm  / sensor_max_x;
        sensor_y_to_mm_factor = sensor_height_mm / sensor_max_y;
        
        break;
      }
      else
      {
        curr_port.stop(); 
      }
    }
    return (serial_port != null);
  }
  
  public void closeConnection()
  {
    setLEDBrightnessAll((byte)0);
    serial_port.stop();
  }
  
  public void setFrameContentControl(byte content)
  {
    senselWriteReg(SENSEL_REG_SCAN_CONTENT_CONTROL, 1, content);
  }
  
  public void setLEDBrightness(int idx, byte brightness)
  {
    if(idx < 16)
      senselWriteReg(SENSEL_REG_LED_BRIGHTNESS + idx, 1, brightness); 
  }
  
  public void setLEDBrightnessAll(byte brightness)
  {
    for(int i = 0; i < 16; i++)
      senselWriteReg(SENSEL_REG_LED_BRIGHTNESS + i, 1, brightness); 
  }
  
  public float getSensorWidthMM()
  {
    return sensor_width_mm;
  }
  
  public float getSensorHeightMM()
  {
    return sensor_height_mm; 
  }
  
  public int getMaxNumContacts()
  {
    return sensor_max_contacts;
  }
  
  public void startScanning()
  {
    senselWriteReg(SENSEL_REG_SCAN_ENABLED, 1, 1);
  }
  
  public void stopScanning()
  {
    senselWriteReg(SENSEL_REG_SCAN_ENABLED, 1, 0);
  }
  
  private int[] _readReg(int addr, int size)
  {
    serial_port.write(SENSEL_READ_HEADER);
    serial_port.write((byte)addr);
    serial_port.write((byte)size);
    
    int[] rx_buf = new int[size]; // TODO (Ilya): I think the rx_buf should be a byte array, and this funciton should return bytes.
    
    int ack;
    while((ack = serial_port.read()) == -1);
    
    if(ack != SENSEL_PT_READ_ACK)
      println("FAILED TO RECEIVE ACK ON READ (regaddr=" + addr + ", ack=" + ack + ")");
    
    int size0;
    while((size0 = serial_port.read()) == -1);
  
    int size1;
    while((size1 = serial_port.read()) == -1);
    
    int resp_size = _convertBytesTo16((byte)size0, (byte)size1);
  
    //if(size != resp_size)
    //  println("RESP_SIZE != SIZE (" + resp_size + "!=" + size + ") ON READ");
    //else
    //  println("RESP_SIZE == SIZE (" + resp_size + "==" + size + ") ON READ");
    
    int checksum = 0;
    
    for(int i = 0; i < size; i++)
    {
       while((rx_buf[i] = serial_port.read()) == -1);
       checksum += rx_buf[i];
    }
    
    checksum = (checksum & 0xFF);
    
    int resp_checksum;
    while((resp_checksum = serial_port.read()) == -1);
    
    if(checksum != resp_checksum)
      println("CHECKSUM FAILED: " + checksum + "!=" + resp_checksum + " ON READ");
    
    return rx_buf;
  }
  
  private int _readContactFrameSize()
  {
    serial_port.write(SENSEL_READ_HEADER);
    serial_port.write((byte)SENSEL_REG_SCAN_READ_FRAME);
    serial_port.write((byte)0x00);
    
    int ack;
    while((ack = serial_port.read()) == -1);
    
    if(ack != SENSEL_PT_FRAME)
      println("FAILED TO RECEIVE FRAME PACKET TYPE ON FRAME READ");
    
    int size0;
    while((size0 = serial_port.read()) == -1);
  
    int size1;
    while((size1 = serial_port.read()) == -1);
    
    int content_bitmask;
    while((content_bitmask = serial_port.read()) == -1);
    //println("CBM: " + content_bitmask);
    
    int frame_counter;
    while((frame_counter = serial_port.read()) == -1);
    //println("FC: " + frame_counter);
    
    //println("Finished reading contact frame size: " + (size0 | (size1 <<8)));  
    
    return _convertBytesTo16((byte)size0, (byte)size1) - 2; //Packet size includes content bitmask and lost frame count which we've already read out
  }
  
  //We only support single-byte writes at this time TODO: Implement multi-byte write
  private void senselWriteReg(int addr, int size, int data)
  {
    if(size != 1)
      println("writeReg only supports writes of size 1");
      
    serial_port.write(SENSEL_WRITE_HEADER);
    serial_port.write((byte)addr);
    serial_port.write((byte)size);
    serial_port.write((byte)data);
    serial_port.write((byte)data); //Checksum
    
    int ack;
    while((ack = serial_port.read()) == -1);
    
    if(ack != SENSEL_PT_WRITE_ACK)
      println("FAILED TO RECEIVE ACK ON WRITE (regaddr=" + addr + ", ack=" + ack + ")");
  }
  
  private int _convertBytesTo32(byte b0, byte b1, byte b2, byte b3)
  {
    return ((((int)b3) & 0xff) << 24) | ((((int)b2) & 0xff) << 16) | ((((int)b1) & 0xff) << 8) | (((int)b0) & 0xff); 
  }
  
  private int _convertBytesTo16(byte b0, byte b1)
  {
    return ((((int)b1) & 0xff) << 8) | (((int)b0) & 0xff); 
  }
  
  // Convert two bytes (which represent a two's complement signed 16 bit integer) into a signed int
  private int _convertBytesToS16(byte b0, byte b1)
  {
    return (((int)b1) << 8) | (((int)b0) & 0xff);
  }
  
  public SenselContact[] readContacts()
  {
    SenselContact[] retval = null;
    int contact_frame_size = _readContactFrameSize() + 1; //For checksum!
  
    //println("CFS: " + contact_frame_size);
  
    if(true)//contact_frame_size > 0)
    {  
      //println("Force frame: " + contact_frame_size);
      byte[] contact_frame = new byte[contact_frame_size];
      
      int aval;
      do
      {
        aval = serial_port.available();
        //println("Aval: " + aval);
        delay(1);
      }
      while(aval < contact_frame_size);
      
      int read_count = serial_port.readBytes(contact_frame);
      
      if(read_count < contact_frame_size)
      {
        println("SOMETHING BAD HAPPENED! (" + read_count + " < " + contact_frame_size + ")");
        exit(); 
      }
      
      int num_contacts = ((int)contact_frame[0]) & 0xff;  
    
      //print("Num Contacts: " + num_contacts + "....");
      
      int idx = 0;
      
      SenselContact[] c = new SenselContact[num_contacts];
      
      for(int i = 0; i < num_contacts; i++)
      {
        c[i] = new SenselContact();
        c[i].total_force = _convertBytesTo32(contact_frame[++idx], contact_frame[++idx], contact_frame[++idx], contact_frame[++idx]);
        c[i].uid = _convertBytesTo32(contact_frame[++idx], contact_frame[++idx], contact_frame[++idx], contact_frame[++idx]);
        //Convert area to square mm
        c[i].area_mm_sq = ((float)_convertBytesTo32(contact_frame[++idx], contact_frame[++idx], contact_frame[++idx], contact_frame[++idx])) * sensor_area_to_mm_sq_factor;
        //Convert x_pos to x_pos_mm
        c[i].x_pos_mm = ((float)_convertBytesTo16(contact_frame[++idx], contact_frame[++idx])) * sensor_x_to_mm_factor;
        //Convert y_pos to y_pos_mm
        c[i].y_pos_mm = ((float)_convertBytesTo16(contact_frame[++idx], contact_frame[++idx])) * sensor_y_to_mm_factor;
        //Convert dx to dx_mm
        c[i].dx_mm =    ((float)_convertBytesTo16(contact_frame[++idx], contact_frame[++idx])) * sensor_x_to_mm_factor;
        //Convert dy to dy_mm
        c[i].dy_mm =    ((float)_convertBytesTo16(contact_frame[++idx], contact_frame[++idx])) * sensor_y_to_mm_factor;
        //Convert orientation to angle in degrees
        c[i].orientation_degrees = ((float)_convertBytesToS16(contact_frame[++idx], contact_frame[++idx])) * sensor_orientation_to_degrees_factor;
        //Convert major_axis to mm (assumes that x_to_mm and y_to_mm are the same)
        c[i].major_axis_mm = ((float)_convertBytesTo16(contact_frame[++idx], contact_frame[++idx])) * sensor_x_to_mm_factor;
        //Convert minor_axis to mm (assumes that x_to_mm and y_to_mm are the same)
        c[i].minor_axis_mm = ((float)_convertBytesTo16(contact_frame[++idx], contact_frame[++idx])) * sensor_x_to_mm_factor;
        ++idx; //peak_x
        ++idx; //peak_y
        c[i].id = (((int)contact_frame[++idx]) & 0xff);
        c[i].type = (((int)contact_frame[++idx]) & 0xff);
      }
      retval = c;
    }  
    //TODO: ACTUALLY USE CHECKSUM!!!
    //byte checksum;
    //while((checksum = (byte)serial_port.read()) == -1);
    //println("finish read");
  
    return retval;
  }
  
  
  // Returns (x,y,z) acceleration in G's using the following coordinate system:
  //
  //          ---------------------------
  //        /   Z /\  _                 /
  //       /       |  /| Y             /
  //      /        | /                /
  //     /         |/                /
  //    /           -----> X        /
  //   /                           /
  //   ----------------------------
  //
  // Assumes accelerometer is configured to the default +/- 2G range
  public float[] readAccelerometerData()
  {
    // Read accelerometer data bytes for X, Y and Z
    int[] acc_bytes = _readReg(SENSEL_REG_ACCEL_X,6);

    // Convert raw bytes to signed values
    int[] acc_values = new int[3];
    acc_values[0] = _convertBytesToS16((byte)acc_bytes[0], (byte)acc_bytes[1]);
    acc_values[1] = _convertBytesToS16((byte)acc_bytes[2], (byte)acc_bytes[3]);
    acc_values[2] = _convertBytesToS16((byte)acc_bytes[4], (byte)acc_bytes[5]);
    
    // Rescale to G's (at a range of +/- 2G, accelerometer returns 0x4000 for 1G acceleration)
    float[] acc_data = new float[3];
    acc_data[0] = ((float)acc_values[0] / 0x4000);
    acc_data[1] = ((float)acc_values[1] / 0x4000);    
    acc_data[2] = ((float)acc_values[2] / 0x4000);
    
    return acc_data;
  }
}
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "Sensel_Morph_Planets" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
