import org.graphper.api.Node;

import java.util.List;

public class NodesDTO {

    private Node node1;

    private List<Node> level1;

    private List<Node> level2;

    public Node getNode1() {
        return node1;
    }

    public void setNode1(Node node1) {
        this.node1 = node1;
    }

    public List<Node> getLevel1() {
        return level1;
    }

    public void setLevel1(List<Node> level1) {
        this.level1 = level1;
    }

    public List<Node> getLevel2() {
        return level2;
    }

    public void setLevel2(List<Node> level2) {
        this.level2 = level2;
    }

    @Override
    public String toString() {
        return "NodesDTO{" +
                "node1=" + node1.toString() +
                ", level1=" + level1.getFirst() +
                ", level2=" + level2.getFirst() +
                '}';
    }
}
