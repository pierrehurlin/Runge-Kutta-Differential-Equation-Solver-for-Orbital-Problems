public class Titan{

	//Maastricht University DKE, 1st Year 
	//Pierre HURLIN, i6253450
	//Titan class, 1st version
	//05/08/2021

	public double x;																//x position
	public double y;																//y position
	public double vx;																//x velocity
	public double vy;																//y velocity

	public Titan(double x0, double y0, double vx0, double vy0){						//Titan constructor

		this.x = x0;
		this.y = y0;
		this.vx = vx0;
		this.vy = vy0;


	}

	public double getx(){															//getx method

		return x;
	}

	public double gety(){															//gety method

		return y;
	}

	public double getvx(){															//getvx method

		return vx;
	}

	public double getvy(){															//getvy method

		return vy;
	}

	public Titan add(Titan other){													//add method

		return new Titan(x+other.getx(), y+other.gety(),vx+other.getvx(),vy+other.getvy());			//adds the x and y positions and velocities of 2 states of Titan
	}

	public Titan mul(double scalar){												//mul method

		return new Titan(x*scalar, y*scalar,vx*scalar,vy*scalar);									//multiply the x and y postions and velocities of a State of Titan by a double 
	}
}