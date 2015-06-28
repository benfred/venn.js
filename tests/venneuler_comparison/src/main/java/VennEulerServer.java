import com.google.gson.Gson;
import com.sun.net.httpserver.Headers;
import edu.uic.ncdm.venn.VennAnalytic;
import edu.uic.ncdm.venn.VennDiagram;
import edu.uic.ncdm.venn.data.VennData;

import java.util.HashMap;
import java.io.IOException;
import java.io.OutputStream;
import java.io.InputStream;
import java.net.InetSocketAddress;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;


public class VennEulerServer {
    class SetSizes {
        public String[] sets;
        public double size;
        SetSizes() { };
    }

    static class Circle {
        public String label;
        public double x;
        public double y;
        public double radius;

        Circle() {}
    }

    public static String jsonFromVennDiagram(VennDiagram vd) {

        Circle[] circles = new Circle[vd.circleLabels.length];
        for (int i = 0; i < circles.length; ++i) {
            circles[i] = new Circle();
            circles[i].radius = vd.diameters[i] / 2;
            circles[i].x = vd.centers[i][0];
            circles[i].y = vd.centers[i][1];
            circles[i].label = vd.circleLabels[i];
        }

        Gson gson = new Gson();
        return gson.toJson(circles);
    }

    public static VennData vennDataFromJSON(String json) {
        Gson gson = new Gson();
        SetSizes[] sizes = gson.fromJson(json, SetSizes[].class);

        String [] regions = new String[sizes.length];
        double [] areas = new double[sizes.length];

        for (int i = 0; i < sizes.length; ++i) {
            StringBuilder s = new StringBuilder();
            for (int j = 0; j < sizes[i].sets.length; ++j) {
                if (j != 0) s.append("&");
                s.append(sizes[i].sets[j]);
            }

            regions[i] = s.toString();
            areas[i] = sizes[i].size;
        }

        VennData data = new VennData(regions, areas);
        return data;
    }

    public static void main(String[] args) throws IOException {
        System.out.println("starting webserver");
        HttpServer server = HttpServer.create(new InetSocketAddress(4242), 0);
        server.createContext("/venneuler", new VennServer());
        server.setExecutor(null);
        server.start();
    }

    static class VennServer implements HttpHandler {
        public void handle(HttpExchange t) throws IOException {
            // read json data from request
            InputStream input = t.getRequestBody();
            java.util.Scanner s = new java.util.Scanner(input).useDelimiter("\\A");
            String serialized =  s.hasNext() ? s.next() : "";

            System.out.printf("Serialized:%s\n", serialized);

            // calculate venn layout
            VennData data = vennDataFromJSON(serialized);
            VennAnalytic va = new VennAnalytic();
            VennDiagram vd = va.compute(data);
            System.out.printf("computed, stress=%.4f\n", vd.stress);
            String output = jsonFromVennDiagram((vd));

            Headers h = t.getResponseHeaders();
            h.set("Access-Control-Allow-Origin", "*");
            h.set("ContentType", "application/json");

            t.sendResponseHeaders(200, output.length());

            OutputStream os = t.getResponseBody();
            os.write(output.getBytes());
            os.close();

        }
    }
}