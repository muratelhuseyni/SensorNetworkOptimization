package SensorNetwork;


public class SensorConnection {
	
	private SensorNode out;
	private SensorNode in;
	private double weight;
	private boolean alive;
	public SensorConnection(SensorNode out, SensorNode in) {
		this.out = out;
		this.in = in;
		this.weight=0.2;
	}

	public boolean isAlive() {
		return alive;
	}

	public void setAlive(boolean alive) {
		this.alive = alive;
	}

	public SensorNode getOut() {
		return out;
	}

	public void setOut(SensorNode out) {
		this.out = out;
	}

	public SensorNode getIn() {
		return in;
	}

	public void setIn(SensorNode in) {
		this.in = in;
	}

	public double getWeight() {
		return weight;
	}

	
	@Override
	public String toString() {
		return "Edge [out=" + out + ", in=" + in + ", weight=" + weight + ", alive=" + alive + "]";
	}

}
