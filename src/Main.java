import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.TreeMap;

import com.jme3.app.SimpleApplication;
import com.jme3.light.AmbientLight;
import com.jme3.light.DirectionalLight;
import com.jme3.material.Material;
import com.jme3.math.ColorRGBA;
import com.jme3.math.Triangle;
import com.jme3.math.Vector3f;
import com.jme3.scene.Geometry;
import com.jme3.scene.Mesh;
import com.jme3.scene.Spatial;
import com.jme3.scene.shape.Sphere;

public class Main extends SimpleApplication{
	
	public static final double THRESHOLD = 1e-6;

	public static void main(String[] args) {
		Main main = new Main();
		main.start();
	}
	
	@Override
	public void simpleInitApp() {
		Spatial mesh = assetManager.loadModel("sphere.obj");
		Material mat = new Material(assetManager, "Common/MatDefs/Light/Lighting.j3md");
        mat.setColor("Diffuse", ColorRGBA.Gray);
        mat.setColor("Ambient", ColorRGBA.White.mult(0.1f));
        mat.setBoolean("UseMaterialColors", true);
        mesh.setMaterial(mat);
        rootNode.attachChild(mesh);
        
        List<Vertex> verticesList = new ArrayList<Vertex>();
        
        if(mesh instanceof Geometry) {
        	Geometry geo = (Geometry) mesh;
        	Mesh m = geo.getMesh();
        	
        	for(int i = 0; i < m.getTriangleCount(); i++) {
        		Triangle t = new Triangle();
        		m.getTriangle(i, t);
        		
        		Vertex v1 = findOrCreateVertex(t.get1(), verticesList);
        		Vertex v2 = findOrCreateVertex(t.get2(), verticesList);
        		Vertex v3 = findOrCreateVertex(t.get3(), verticesList);
        		
        		if(!v1.neighbors.contains(v2)) v1.neighbors.add(v2);
        		if(!v1.neighbors.contains(v3)) v1.neighbors.add(v3);
        		if(!v2.neighbors.contains(v1)) v2.neighbors.add(v1);
        		if(!v2.neighbors.contains(v3)) v2.neighbors.add(v3);
        		if(!v3.neighbors.contains(v1)) v3.neighbors.add(v1);
        		if(!v3.neighbors.contains(v2)) v3.neighbors.add(v2);
        	}
        	
        	Vertex seed = new Vertex();
        	seed.pos = new Vector3f(0f, -1f, 0f);
        	computeGeodesicDistance(verticesList, seed);
        	
        	TreeMap<Double, List<Vertex>> rows = groupByRows(verticesList);
        	List<Vertex> r = rows.get(rows.higherKey(1d));
        	
        	for(Vertex v : r) {
        		Material marker = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
            	marker.setColor("Color", ColorRGBA.Red);
            	Geometry mg = new Geometry("VertexMarker", new Sphere(8, 8, 0.02f));
            	mg.setMaterial(marker);
            	mg.setLocalTranslation(v.pos);
            	rootNode.attachChild(mg);
        	}
        	
        	System.out.println(verticesList.size());
        	
//        	VertexBuffer vb = m.getBuffer(VertexBuffer.Type.Position);
//        	float[] vertices = BufferUtils.getFloatArray((FloatBuffer) vb.getData());
//        	
//        	Material marker = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
//        	marker.setColor("Color", ColorRGBA.Red);
//        	
//        	// For each vertex, create a small sphere and position it at the vertex's location
//            for (int i = 0; i < vertices.length; i += 3) {
//                Vector3f vertexPos = new Vector3f(vertices[i], vertices[i + 1], vertices[i + 2]);
//                Geometry mg = new Geometry("VertexMarker", new Sphere(8, 8, 0.01f)); // Small sphere
//                mg.setMaterial(marker);
//                mg.setLocalTranslation(vertexPos);
//                rootNode.attachChild(mg);
//            }
        }

        // Set up lighting
        AmbientLight ambient = new AmbientLight();
        ambient.setColor(ColorRGBA.White);
        rootNode.addLight(ambient);

        DirectionalLight sun = new DirectionalLight();
        sun.setDirection(new Vector3f(-1, -1, -1).normalizeLocal());
        sun.setColor(ColorRGBA.White);
        rootNode.addLight(sun);
        
        DirectionalLight backLight = new DirectionalLight();
        backLight.setDirection(new Vector3f(1, 1, 1).normalizeLocal());
        backLight.setColor(ColorRGBA.White);
        rootNode.addLight(backLight);

        cam.setLocation(new Vector3f(0, 0, 10));
        cam.lookAt(Vector3f.ZERO, Vector3f.UNIT_Y);
	}
	
	public Vertex findOrCreateVertex(Vector3f pos, List<Vertex> verticesList) {
		for(Vertex v : verticesList) {
			if(v.pos.equals(pos)) {
				return v;
			}
		}
		
		Vertex newVertex = new Vertex();
		newVertex.pos = pos;
		verticesList.add(newVertex);
		return newVertex;
	}
	
	public void computeGeodesicDistance(List<Vertex> vertices, Vertex seed) {
		PriorityQueue<Vertex> queue = new PriorityQueue<>(Comparator.comparing(v -> v.f));
		
		seed.f = 0;
		queue.add(seed);
		
		while(!queue.isEmpty()) {
			Vertex current = queue.poll();
			current.finalized = true;
			
			for(Vertex neighbor : current.neighbors) {
				if(!neighbor.finalized) {
					double newDist = current.f + current.pos.distance(neighbor.pos);
					
					if(newDist < neighbor.f) {
						neighbor.f = newDist;
						queue.remove(neighbor);
						queue.add(neighbor);
					}
				}
			}
		}
	}
	
	public static boolean areClose(double a, double b) {
		return Math.abs(a - b) < THRESHOLD;
	}
	
	public static TreeMap<Double, List<Vertex>> groupByRows(List<Vertex> vertices) {
		TreeMap<Double, List<Vertex>> rows = new TreeMap<>();
		
		for(Vertex vertex : vertices) {
			boolean added = false;
			
			for(Double key : rows.keySet()) {
				if(areClose(vertex.f, key)) {
					rows.get(key).add(vertex);
					added = true;
					break;
				}
			}
			
			if(!added) {
				List<Vertex> newRow = new ArrayList<>();
				newRow.add(vertex);
				rows.put(vertex.f, newRow);
			}
		}
		
		return rows;
	}
}