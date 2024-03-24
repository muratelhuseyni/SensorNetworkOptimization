package SensorNetwork;

import java.util.ArrayList;


public class SensorNode implements Comparable<SensorNode>{
	private String id;
	private String address;
	private int port;
	private int load;
	private int infectionscore;
	private SensorInfo sensorInfo;
	private SensorType type;
	private ArrayList<SensorConnection> sensorConnections;
	private ArrayList<SensorNode> neighbours;
	private ArrayList<SensorConnection> inlist;
	private ArrayList<SensorConnection> outlist;
	private double fv;
	
	
	public SensorNode(String id, String type, String address, int port) {
		this.id=id;
		this.address=address;
		this.port=port;
		this.load=0;
		this.sensorConnections =new ArrayList<SensorConnection>();
		this.neighbours=new ArrayList<SensorNode>();
		this.inlist = new ArrayList<SensorConnection>();
		this.outlist= new ArrayList<SensorConnection>();
	}
	
	public SensorNode(String id) {
		this.id=id;
		this.sensorConnections =new ArrayList<SensorConnection>();
		this.neighbours=new ArrayList<SensorNode>();
		this.inlist = new ArrayList<SensorConnection>();
		this.outlist= new ArrayList<SensorConnection>();
		this.load=0;
	}
		
	public SensorType getNodeFormat() {
		return type;
	}

	public void setType(String type) {
		if(type.equals("Sensor")) {
			this.type= SensorType.SENSOR;
		}else if(type.equals("Gateway")) {
			this.type= SensorType.GATEWAY;
		}else {
			System.out.println("SensorID:"+this.id+"Hib�s szenzort�pus, a program kil�p.");
			System.exit(0);
		}
	}
	
	public SensorNode getminpackageneighbour() {
		SensorNode min=this.neighbours.get(0);
		for(int i=1;i<this.neighbours.size();i++) {
			if(this.neighbours.get(i).getLoad()<min.getLoad()) {
				min=this.neighbours.get(i);
			}
		}
		return min;
	}
	
	public int getLoad() {
		return load;
	}

	public void addPackagetoNode() {
		this.load++;
	}
	
	public void setLoad(int n) {
		this.load=n;
	}

	public void setType(SensorType type) {
		this.type = type;
	}

	public void setEdges(ArrayList<SensorConnection> sensorConnections) {
		this.sensorConnections = sensorConnections;
	}

	public String getId() {
		return id;
	}
	
	public void setId(String id) {
		this.id = id;
	}
	
	public String getAddress() {
		return address;
	}
	
	public void setAddress(String address) {
		this.address = address;
	}
	
	public int getPort() {
		return port;
	}
	
	public void setPort(int port) {
		this.port = port;
	}
	
	public SensorInfo getState() {
		return sensorInfo;
	}
	
	public void setState(SensorInfo sensorInfo) {
		this.sensorInfo = sensorInfo;
	}
	
	public ArrayList<SensorNode> getNeighbours() {
		return neighbours;
	}
	
	public void setNeighbours(ArrayList<SensorNode> neighbours) {
		this.neighbours = neighbours;
	}
	
	public void addNeighbour(SensorNode a) {
		this.neighbours.add(a);
	}
	
	public void addOutEdge(SensorConnection e) {
		this.outlist.add(e);
	}
	
	public void addinEdge(SensorConnection e) {
		this.inlist.add(e);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((address == null) ? 0 : address.hashCode());
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		result = prime * result + port;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SensorNode other = (SensorNode) obj;
		if (address == null) {
			if (other.address != null)
				return false;
		} else if (!address.equals(other.address))
			return false;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		if (port != other.port)
			return false;
		return true;
	}

	@Override
	public int compareTo(SensorNode n) {
		return this.load > n.getLoad() ? -1 : (this.load < n.getLoad() ?  1 : 0);
	}	
	
	public void finalizefv(int samplesize) {
		this.fv=this.fv/(double)samplesize;
	}
	

	@Override
	public String toString() {
		return "Node [id=" + id + ", address=" + address + ", port=" + port + "]";
	}

	public ArrayList<SensorConnection> getEdges() {
		return sensorConnections;
	}

	public void addEdge(SensorConnection sensorConnection) {
		this.sensorConnections.add(sensorConnection);
	}

	public double getFv() {
		return fv;
	}

	public void setFv(double fv) {
		this.fv = fv;
	}
	
	public void addFv(int i) {
		this.fv+=i;
	}

	public ArrayList<SensorConnection> getInlist() {
		return inlist;
	}

	public void setInlist(ArrayList<SensorConnection> inlist) {
		this.inlist = inlist;
	}

	public ArrayList<SensorConnection> getOutlist() {
		return outlist;
	}

	public void setOutlist(ArrayList<SensorConnection> outlist) {
		this.outlist = outlist;
	}

	public int getInfectionscore() {
		return infectionscore;
	}

	public void setInfectionscore(int infectionscore) {
		this.infectionscore = infectionscore;
	}

}
