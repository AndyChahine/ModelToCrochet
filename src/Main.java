import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
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
	
	public static final double THRESHOLD = 1e-3;
	public static final double GRID_WIDTH = 0.15d;

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
        	
        	computeGeodesicDistance(verticesList, verticesList.get(0));
        	
        	Vertex ver = new Vertex();
        	ver.pos = new Vector3f(0, 0, 0);
        	Vertex v = findClosestVertex(ver, verticesList);
        	
        	TreeMap<Double, List<Vertex>> rows = groupByRows(verticesList);
        	resampleVertices(rows);
        	
        	for(Map.Entry<Double, List<Vertex>> entry : rows.entrySet()) {
        		List<Vertex> r = entry.getValue();
        		System.out.println(r.size() + " " + r);
        		
        		for(int i = 0; i < r.size(); i++) {
            		Material marker = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
                	marker.setColor("Color", new ColorRGBA(1f, (float)i / (float)r.size(), (float)i / (float)r.size(), 1f));
                	Geometry mg = new Geometry("VertexMarker", new Sphere(8, 8, 0.02f));
                	mg.setMaterial(marker);
                	mg.setLocalTranslation(r.get(i).pos);
                	rootNode.attachChild(mg);
            	}
        	}
        	
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
        flyCam.setMoveSpeed(10f);
        flyCam.setZoomSpeed(10f);
        flyCam.setDragToRotate(true);
        inputManager.setCursorVisible(true);
	}
	
	public void resampleVertices(TreeMap<Double, List<Vertex>> treeMap) {
		double epsilon = GRID_WIDTH * 0.0001d;
		
		for(Entry<Double, List<Vertex>> rowEntry : treeMap.entrySet()) {
			
			Double key = rowEntry.getKey();
			List<Vertex> rowVertices = rowEntry.getValue();
			
			List<Vertex> selectedVertices = new ArrayList<>();
			selectedVertices.add(rowVertices.get(0));
			
			for(int i = 1; i < rowVertices.size(); i++) {
				double dist = rowVertices.get(i).pos.distance(selectedVertices.get(selectedVertices.size() - 1).pos);
				
				if(Math.abs(dist - GRID_WIDTH) <= epsilon) {
					selectedVertices.add(rowVertices.get(i));
				}else if(dist > GRID_WIDTH) {
					Vertex interVertex = interpolate(selectedVertices.get(selectedVertices.size() - 1), rowVertices.get(i), GRID_WIDTH);
					selectedVertices.add(interVertex);
					i--;
				}
			}
			
			treeMap.put(key, selectedVertices);
		}
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
	
	public static Vertex findClosestVertex(Vertex input, List<Vertex> vertices) {
		
		if(vertices == null || vertices.isEmpty()) {
			return null;
		}
		
		Vertex closest = vertices.get(0);
		double minDist = input.pos.distance(closest.pos);
		for(Vertex vertex : vertices) {
			
			double curDist = input.pos.distance(vertex.pos);
			if(curDist < minDist) {
				minDist = curDist;
				closest = vertex;
			}
		}
		
		return closest;
	}
	
	public static Vertex interpolate(Vertex start, Vertex end, double targetDist) {
		Vector3f direction = end.pos.subtract(start.pos);
		
		direction.normalizeLocal();
		
		direction.multLocal((float) targetDist);
		
		Vector3f interpolatedPos = start.pos.add(direction);
		
		Vertex interpolatedVertex = new Vertex();
		interpolatedVertex.pos = interpolatedPos;
		
		return interpolatedVertex;
	}
}