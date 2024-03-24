package DataProvider;

import java.io.FileReader;
import java.io.IOException;

import com.opencsv.CSVReader;

import SensorNetwork.*;

public class DataReader {
    public static SensorNetwork readEdges(String filepath) throws IOException {
        SensorNetwork sensorNetwork = new SensorNetwork();

        try (CSVReader reader = new CSVReader(new FileReader(filepath), ';')) {
            String[] line;

            reader.readNext();

            while ((line = reader.readNext()) != null) {

                if (!line[0].equals(line[1])) {
                    sensorNetwork.addEdge(new SensorConnection(new SensorNode(line[0]), new SensorNode(line[1])));
                }
            }
        }

        return sensorNetwork;
    }

    public static SensorNetwork readNodes(String filepath, SensorNetwork sensorNetwork) throws IOException {
        int processedNodes = 0;
        try (CSVReader reader = new CSVReader(new FileReader(filepath), ';')) {
            reader.readNext();

            String[] line;
            while ((line = reader.readNext()) != null) {
                processedNodes++;
                String nodeId = line[0];
                SensorNode sensorNode = sensorNetwork.getNodeMap().get(nodeId);

                if (sensorNode != null) {
                    sensorNode.setType(line[1]);
                    sensorNode.setAddress(line[2]);
                    sensorNode.setPort(Integer.parseInt(line[3]));

                    switch (sensorNode.getNodeFormat()) {
                        case SENSOR:
                            sensorNetwork.addSensor(sensorNode);
                            break;
                        case GATEWAY:
                            sensorNetwork.addGateway(sensorNode);
                            break;
                    }
                }
            }
        }

        if (processedNodes != sensorNetwork.getNodes().size()) {
            System.err.println("Mismatch in the number of nodes processed vs. expected.");
        }

        return sensorNetwork;
    }
}
