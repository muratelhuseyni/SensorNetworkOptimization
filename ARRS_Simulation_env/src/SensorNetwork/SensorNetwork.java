package SensorNetwork;

import java.util.*;


public class SensorNetwork {
	private Map<String, SensorNode> nodemap; // Map of nodes to search by index;
	private Set<SensorNode> sensorNodes; //Set of nodes
	private ArrayList<SensorConnection> sensorConnections; //List of edges
	private ArrayList<SensorNode> sensors;
	private ArrayList<SensorNode> gateways;

	
	
	public SensorNetwork() {
		this.nodemap=new HashMap<String, SensorNode>();
		this.sensorNodes =new LinkedHashSet<SensorNode>();
		this.sensorConnections =new ArrayList<SensorConnection>();
		this.sensors=new ArrayList<SensorNode>();
		this.gateways=new ArrayList<SensorNode>();
	}
	
	
	
	public ArrayList<SensorNode> getSensors() {
		return sensors;
	}

	

	public void addSensor(SensorNode sensor) {
		this.sensors.add(sensor);
	}



	public ArrayList<SensorNode> getGateways() {
		return gateways;
	}

	public void removeGatewayFromSensorList(SensorNode n){
		this.sensors.remove(n);
	}

	public void removeSensorFromGatewayList(SensorNode n){
		this.gateways.remove(n);
	}

	public void addGateway(SensorNode gateway) {
		this.gateways.add(gateway);
	}



	public void Addnode(SensorNode sensorNode) {
		sensorNodes.add(sensorNode);
	}
	
	public Map<String, SensorNode> getNodeMap(){
		return nodemap;
	}

	public void addNode(SensorNode n) {
		sensorNodes.add(n);
		nodemap.put(n.getId(), n);
	}
	
	public Set<SensorNode> getNodes(){
		return sensorNodes;
	}
	
	public void addEdge(SensorConnection a) {
		SensorNode in=null;
		SensorNode out=null;
	
		
		if(sensorNodes.contains(a.getOut())) {
			out=(SensorNode)nodemap.get(a.getOut().getId());
		}else {
			sensorNodes.add(a.getOut());
			out=a.getOut();
		}
		
		if(sensorNodes.contains(a.getIn())) {
			in=(SensorNode)nodemap.get(a.getIn().getId());
		}else {
			sensorNodes.add(a.getIn());
			in=a.getIn();
		}
		a.setIn(in);
		a.setOut(out);
		//UNDIRECTED
		out.addNeighbour(in);
		in.addNeighbour(out);
		
		out.addOutEdge(a);
		in.addinEdge(a);
		in.addEdge(a);
		out.addEdge(a);
		sensorConnections.add(a);
		nodemap.put(in.getId(), in);
		nodemap.put(out.getId(),out);
	}
	
	public ArrayList<SensorConnection> getEdges(){
		return sensorConnections;
	}

	public Map<String, SensorNode> getNodemap() {
		return nodemap;
	}

	public void setNodemap(Map<String, SensorNode> nodemap) {
		this.nodemap = nodemap;
	}

	public void setNodes(Set<SensorNode> sensorNodes) {
		this.sensorNodes = sensorNodes;
	}

	public void setEdges(ArrayList<SensorConnection> sensorConnections) {
		this.sensorConnections = sensorConnections;
	}

	public void empty() {
		for(SensorNode n: sensorNodes) {
			n.setState(SensorInfo.NOTVISITED);
		}
	}

	public void deleteFv() {
		for(SensorNode n: sensorNodes) {
			n.setFv(0);
		}
	}
	
	public double getFvSum() {
		double sum=0.0;
		for(SensorNode n: sensorNodes) {
			sum+=n.getFv();
		}
		return sum;
	}
	
	public void setAllNodeToSensor() {
		this.sensors=new ArrayList<SensorNode>(this.sensorNodes);
		for(SensorNode n: this.gateways) {
			n.setType(SensorType.SENSOR);
		}
		this.gateways= new ArrayList<SensorNode>();
	}

	public long getLoadSum() {
		long sum=0;
		for(SensorNode n: this.getNodes()) {
			sum+=n.getLoad();
		}
		return sum;
	}

	public long getMaxLoad() {
		long maxload=0;
		for(SensorNode n: this.getNodes()) {
			if(maxload<n.getLoad()) {
				maxload=n.getLoad();
			}
		}
		return maxload;
	}
	
	public long getBalance() {
		long average=0;
		for(SensorNode n: this.getNodes()) {
			average+=n.getLoad();
		}
		average=average/this.getNodes().size();
		long objvalue=0;
		for(SensorNode n:this.getNodes()) {
			objvalue=objvalue + Math.abs(average-n.getLoad());
		}
		return objvalue;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((sensorConnections == null) ? 0 : sensorConnections.hashCode());
		result = prime * result + ((nodemap == null) ? 0 : nodemap.hashCode());
		result = prime * result + ((sensorNodes == null) ? 0 : sensorNodes.hashCode());
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
		SensorNetwork other = (SensorNetwork) obj;
		if (sensorConnections == null) {
			if (other.sensorConnections != null)
				return false;
		} else if (!sensorConnections.equals(other.sensorConnections))
			return false;
		if (nodemap == null) {
			if (other.nodemap != null)
				return false;
		} else if (!nodemap.equals(other.nodemap))
			return false;
		if (sensorNodes == null) {
			if (other.sensorNodes != null)
				return false;
		} else if (!sensorNodes.equals(other.sensorNodes))
			return false;
		return true;
	}
	
	@Override
	public String toString() {
		return "Network [nodes=" + sensorNodes.size() + ", edges=" + sensorConnections.size() + "]";
	}

}
