import java.lang.*;
import javax.swing.*;
import java.awt.geom.*;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Color;

public class Compo extends JComponent{

  //Maastricht University DKE, 1st Year 
  //Pierre HURLIN, i6253450
  //Compo class, 1st version
  //05/08/2021

  private static Titan[] z;                                                   //stores the experimental states of Titan
  private static double[] x_true;                                             //stores the theoritical x positions of Titan
  private static double[] y_true;                                             //stores the theoritical y positions of Titan
  private Color color = Color.BLACK;

    public void setColor (Color color) {                                      //setColor method
      this.color = color;
  }

  public void setTitan (Titan[] z) {                                          //setTitan method
    this.z = z;
  }
  public void setTrue (double[] x_true,double[] y_true) {                     //setTrue method
    this.x_true = x_true;
    this.y_true = y_true;
  }

  public void paintComponent(Graphics g) {                                    //paintComponent method

      Graphics2D g2=(Graphics2D)g;

      g2.setColor(Color.BLUE);
      Ellipse2D.Double pos1 = new Ellipse2D.Double(200,200,10,10);            //Focus S
      g2.draw(pos1);

      for (int i = 0;i<z.length;i++){                                         //for loop to draw all the positions of Titan

        g2.setColor(Color.BLACK);
        Ellipse2D.Double pos = new Ellipse2D.Double((z[i].getx()*100)+200,(z[i].gety()*100)+200,5,5);//3rd order Runge-Kutta method
        g2.draw(pos);
        g2.setColor(Color.RED);
        Ellipse2D.Double pos2 = new Ellipse2D.Double((x_true[i]*100)+200,(y_true[i]*100)+200,5,5);//exact solution
        g2.draw(pos2);
      }

  

  }


}
