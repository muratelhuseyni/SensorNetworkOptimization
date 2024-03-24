package SimulationLevel;

import SensorNetwork.SensorNode;
import SensorNetwork.SensorType;

import java.util.ArrayList;
import java.util.List;

public class Packages {
    public static class WRWPackage extends PackageBase {

        public WRWPackage(SensorNode origin) {
            super(origin);
            // TODO Auto-generated constructor stub
        }

        @Override
        public void PackageSim() {
            double weightsum=0.0;
            for(SensorNode n: actnode.getNeighbours()) {
                weightsum+= n.getInfectionscore();
            }
            double random= Math.random()*weightsum;
            double countweight = 0.0;
            for(SensorNode n: actnode.getNeighbours()) {
                countweight += n.getInfectionscore();
                if(countweight>=random) {
                    this.actnode=n;
                    this.actnode.addPackagetoNode();
                    return;
                }
            }

        }

    }

    public static class RWPackage extends PackageBase {

        public RWPackage(SensorNode actnode) {
            super(actnode);
        }

        @Override
        public void PackageSim() {
            this.actnode = this.actnode.getNeighbours().get((int)(Math.random()*(this.actnode.getNeighbours().size()-1)));
            this.actnode.addPackagetoNode();
        }

    }

    public static class MLRWPackage extends PackageBase {

        public MLRWPackage(SensorNode origin) {
            super(origin);
        }

        @Override
        public void PackageSim() {
            this.actnode=this.actnode.getminpackageneighbour();
            this.actnode.addPackagetoNode();
        }

    }

    public static class LWRWPackage extends PackageBase {

        public LWRWPackage(SensorNode origin) {
            super(origin);
        }

        @Override
        public void PackageSim() {
            List<SensorNode> neighbours = actnode.getNeighbours();

            List<Double> weights = new ArrayList<>();
            for (SensorNode n : neighbours) {
                weights.add(1.0 / (n.getLoad() + 1));
            }

            List<Double> packagessum = new ArrayList<>();
            double sum = 0;
            for (double weight : weights) {
                sum += weight;
                packagessum.add(sum);
            }

            double random = Math.random() * sum;
            for (int i = 0; i < neighbours.size(); i++) {
                if (packagessum.get(i) > random) {
                    this.actnode = neighbours.get(i);
                    if (this.actnode.getNodeFormat() == SensorType.SENSOR) { // Gateways do not get any load
                        this.actnode.addPackagetoNode();
                    }
                    return;
                }
            }
        }
    }
}
