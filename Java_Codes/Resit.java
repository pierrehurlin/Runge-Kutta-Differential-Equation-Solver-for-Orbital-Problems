import java.lang.*;
import javax.swing.*;
import java.awt.geom.*;
import java.awt.Graphics2D;
import java.awt.BorderLayout;
import java.awt.Color;

public class Resit{

	//Maastricht University DKE, 1st Year 
	//Pierre HURLIN, i6253450
	//Resit class, 1st version
	//05/08/2021


	private static double e = 0.5;											//Eccentricity
	private static double mu = 1.0;											//Product of the gravitational constant by the mass of Saturn
	
	private static double a = 1.0;											//Semi-major axis
	private static double T = 2*Math.PI*Math.pow((Math.pow(a,3))/mu,1/2);	//Period
	private static double b =  a*Math.pow(1-Math.pow(e,2), 1/2);			//Semi-minor axis
	private static double c = Math.pow(Math.pow(a,2)-Math.pow(b,2), 1/2);	//Distance from the center to the focus S
	private static double q = a*(1-e);										//Distance from the focus to the perikorn
	
	private static double t0 = 0;											//Initial Time
	private static double tf = T;											//Final Time
	private static double h = 0.1;											//Step size
	private static int n = (int)(Math.floor((tf-t0)/h)+1);					//number of steps
	
	private static double c2 = 1.0/3.0;										//Parameters for Heun's method (change for any other 3rd RK method)
	private static double c3 = 2.0/3.0;

	private static double b3 = (3*c2-2)/(6*c3*(c2-c3));						//Constrained parameters for 3rd RK method
	private static double b2 = (3*c3-2)/(6*c2*(c3-c2));
	private static double b1 = 1-b2-b3;
	private static double a21 = c2;
	private static double a32 = 1/(6*c2*b3);
	private static double a31 = c3-a32;

	public Titan newtonlaw(double t, Titan z){								//newton's law of universal gravitation f(t,z)

		double[] f = new double[4];

		f[0] = z.getvx();
		f[1] = z.getvy();
		f[2] = -mu/(Math.pow(Math.pow(z.getx(),2)+Math.pow(z.gety(),2),3.0/2.0))*z.getx();
		f[3] = -mu/(Math.pow(Math.pow(z.getx(),2)+Math.pow(z.gety(),2),3.0/2.0))*z.gety();


		return new Titan(f[0],f[1],f[2],f[3]);
	}


	public static void main(String[] args) {

		Resit resit = new Resit();
	
		double x0 = q;																			//Initial conditions
		double y0 = 0;
		double vx0 = 0;
		double vy0 = Math.pow(mu/a*(1+e)/(1-e),1.0/2.0);

		Titan z = new Titan(x0,y0,vx0,vy0);
		Titan[] res = new Titan[n];
		res = resit.RKsolve(n,z);

		double[] u = resit.error(res,0);														//eccentric anomaly
		double[] x_true = resit.xtrue(u);														//Theoritical x position of Titan at each time-step
		double[] y_true = resit.ytrue(u);														//Theoritical y position of Titan at each time-step

		double[] re_x = resit.rex(res,x_true);													//relative error on x at each time-step
		double[] re_y = resit.rey(res,y_true);													//relative error on x at each time-step
		double rextot = 0;
		double reytot = 0;

		for (int j = 0;j<u.length;j++){

			System.out.println("i = " + j + " x = "+ res[j].getx() +" y = " + res[j].gety() + " vx = " + res[j].getvx() + " vy = " + res[j].getvy());		//prints the state of Titan at each time step
			System.out.println("xtrue = " + x_true[j] + " ytrue = " + y_true[j]);																			//prints the theoretical position of Titan
			/*System.out.println("err = " + u[j]);																											//eccentric anomaly
			System.out.println("rex = " + re_x[j]);																											//prints the relative error on x
			System.out.println("rey = " + re_y[j]);*/																										//prints the relative error on y

			if (y_true[j] != 0.0){
				reytot = reytot + re_y[j];														//sum of the relative error on x
			}	
			if (x_true[j]!=0.0){
				rextot = rextot + re_x[j];														//sum of the relative error on y
			}




		}
		double rexave = rextot/re_x.length;														//average relative error on x
		double reyave = reytot/re_y.length;														//average relative error on y
		System.out.println("step size: " + h);
		System.out.println(rexave);
		System.out.println(reyave);

		JFrame frame = new JFrame();
		frame.setLayout(new BorderLayout());													
		Compo component = new Compo();
		component.setTitan(res);																//GUI
		component.setTrue(x_true,y_true);
		frame.add(component,BorderLayout.CENTER);
		frame.setTitle("Graph");
		frame.setSize(400,400);
		frame.setLocationRelativeTo(null);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    	frame.setVisible(true);




	}

	public Titan[] RKsolve(int n, Titan z){													//RKsolve method

		Titan[] res = new Titan[n+1];
		res[0] = z;
		double t = t0;
		//Titan j = newtonlaw(t,titan);

		//System.out.println("x = "+ j.getx() +" y = " + j.gety() + " vx = " + j.getvx() + " vy = " + j.getvy());

		for (int i = 1; i<n+1;i++){															//iterates the Rkstep method for n times

			t = t+h*i;

			res[i] = RKstep(res[i-1],t);
		}

		return res;																			//returns an array composed of each states of Titan during tf time

	}

	public Titan RKstep(Titan z,double t){													//RK 3rd order step method

		Titan k1 = newtonlaw(t, z).mul(h);
		Titan k2 = newtonlaw(t+c2*h,z.add(k1.mul(a21))).mul(h);
		Titan k3 = newtonlaw(t+c3*h,z.add(k1.mul(a31)).add(k2.mul(a32))).mul(h);
		Titan w_RK = z.add(k1.mul(b1)).add(k2.mul(b2)).add(k3.mul(b3));

		return w_RK;																		//returns new state of Titan after the time step


	}

	public double[] error(Titan[] res,double t0){											//error method 

		double[] u = new double[res.length];												//array that stores the eccentric anomalies

		for (int i = 0;i<u.length;i++){

			double t= t0+h*i;																//time variable
			double Ma = 2*Math.PI/T*t;														//mean anomaly
			double u0 = Ma/(1-e);															//initialisation of uO
			double diff = 1.0;																//Resolution of the Kepler's equation

			while(Math.abs(diff)>0.0001){													//Precision criteria for the resolution

				double u1 = e*Math.sin(u0)+Ma;												//Recurrence equation
				diff=u1-u0;																	//Difference between u1 and u0
				u0 = u1;																	//Previous value
			}
			u[i] = u0;																		//store the eccentric anomaly at this time step
		}
		return u;																			//returns the eccentric anomalies for each time step 
	}

	public double[] xtrue(double[] u){														//xtrue method

		double[] x_true = new double[u.length];

		for (int i = 0;i<u.length;i++){
			x_true[i] = a*(Math.cos(u[i])-e);												//calculates the true x position of the time step
		}

		return x_true;																		//returns the true x positions for each time step

	}
	public double[] ytrue(double[] u){														//ytrue method

		double[] y_true = new double[u.length];

		for (int i = 0;i<u.length;i++){
			y_true[i] = a*Math.sin(u[i])*Math.pow(1-Math.pow(e,2),0.5);						//calculates the true y position of the time step
		}

		return y_true;																		//returns the true y positions for each time step

	}
	public double[] rex(Titan[] res,double[] x_true){										//rex method

		double[] re_x = new double[res.length];

		for (int i = 0;i<res.length;i++){												
			/*System.out.println("getx " + res[i].getx());
			System.out.println("xtrue " + x_true[i]);*/
			re_x[i] = 100*Math.abs(res[i].getx()-x_true[i])/Math.abs(x_true[i]);			//relative error on x for this time step
		}

		return re_x;																		//returns the relative error on x for each time step

	}
	public double[] rey(Titan[] res,double[] y_true){										//rey method

		double[] re_y = new double[res.length];

		for (int i = 0;i<res.length;i++){
			re_y[i] = 100*Math.abs(res[i].gety()-y_true[i])/Math.abs(y_true[i]);			//relative error on y for this time step
		}

		return re_y;																		//returns the relative error on y for each time step


	}




}