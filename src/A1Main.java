
import java.io.PrintWriter;
import java.io.FileWriter;
import java.util.Map;
import java.util.TreeMap;
import java.awt.FlowLayout;
import java.awt.Container;
import java.awt.Button;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Dimension;
import java.awt.Color;
import java.awt.Polygon;
import ch.reto_hoehener.japng.Apng;
import ch.reto_hoehener.japng.ApngFactory;
import ch.reto_hoehener.japng.JapngException;
import java.awt.Label;
import java.util.Collections;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;
import static java.lang.Math.log10;

// Provides a source of random numbers for mutations and crossover
class Rand {

    private static Random random = new Random();

    public static int nextInt(int n) {
        return random.nextInt(n);
    }

    public static float nextFloat() {
        return random.nextFloat();
    }
}

// Represents the x,y coordinates that determine a shape
class Point {

    private int x, y;

    public Point(int x, int y) {
        this.x = x;
        this.y = y;
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }
}

// Represents the geometry of a polygon
// This is the interface for various shapes
// You can add classes that represent other shapes
class Poly {

    private Polygon polygon;
    private List<Point> points;
    private int intWidth, intHeight;

    public Poly(int width, int height) {
        intHeight = height;
        intWidth = width;
        points = new ArrayList<Point>();
        // This code creates a polygon with points chosen uniformly randomly
        // from the area of the original image.
        // You can change this to produce a different distribution of polygons
        int randNum;
        int randDistRate = 100;

        //3 point Polygon
        for (int i = 0; i < 3; i++) {

            randNum = new Random().nextInt(randDistRate);

            if (randNum < randDistRate / 6) {
                points.add(new Point(Rand.nextInt(width / 2), Rand.nextInt(height)));
            } else if (randNum > randDistRate / 6 && randNum < randDistRate / 3) {
                points.add(new Point(Rand.nextInt(width), Rand.nextInt(height / 2)));
            } else if (randNum > randDistRate / 3 && randNum < randDistRate / 2) {
                points.add(new Point(Rand.nextInt(width / 2), Rand.nextInt(height / 2)));
            } else {
                points.add(new Point(Rand.nextInt(width), Rand.nextInt(height)));
            }
        }
        getPolygon();
    }

    public Poly(List<Point> points) {
        this.points = points;
        getPolygon();
    }

    public Poly(Poly p) {
        p.getPolygon();
        polygon = new Polygon(p.polygon.xpoints, p.polygon.ypoints, p.polygon.npoints);
        points = new ArrayList<Point>(p.points);
        intWidth = p.intWidth;
        intHeight = p.intHeight;
    }

    public Polygon getPolygon() {
        if (polygon != null) {
            return polygon;
        }
        polygon = new Polygon();
        for (Point point : points) {
            polygon.addPoint(point.getX(), point.getY());
        }
        return polygon;
    }

    public Polygon getPolygon(int scale) {
        Polygon polygon = new Polygon();
        for (Point point : points) {
            polygon.addPoint(scale * point.getX(), scale * point.getY());
        }
        return polygon;
    }

    void mutate() {

        Polygon p = getPolygon();
        final int MaxPoints = 8;
        final int minPoints = 3;
        final int addPointRate = 3000;
        final int movePointRate = 1000;
        final int mutatePointRate = 2000;

        //Add Point
        if (Rand.nextInt(addPointRate) == 1) {
            if (p.npoints < MaxPoints) {

                p.addPoint(Rand.nextInt(intWidth), Rand.nextInt(intHeight));
            }
        }

        //Move Polygon
        if (Rand.nextInt(movePointRate) == 1) {
            p.translate(Rand.nextInt(intWidth), Rand.nextInt(intHeight));
        }

        //Mutate Point
        if (Rand.nextInt(mutatePointRate) == 1) {

            int randDistRate = new Random().nextInt(100);
            int randNum = new Random().nextInt(randDistRate);


            if (randNum < randDistRate / 6) {
                p.xpoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intWidth / 2);
                p.ypoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intHeight / 2);
            } else if (randNum > randDistRate / 6 && randNum < randDistRate / 3) {
                p.xpoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intWidth / 2);
                p.ypoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intHeight);
            } else if (randNum > randDistRate / 3 && randNum < randDistRate / 2) {
                p.xpoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intWidth);
                p.ypoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intHeight / 2);
            } else {
                p.xpoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intWidth);
                p.ypoints[Rand.nextInt(p.npoints)] = Rand.nextInt(intHeight);
            }
        }
    }
}

// This is the interface for various shapes
// You can add classes that represent other shapes
interface Shape {

    public void draw(Graphics g, int scale);

    public Shape mutate();
}

// Represents a polygon with a specific geometry and color
// This is the interface to a factory that creates shapes
// You can add another shape factory if you add a class that represents a
// new type of shape
// Factory class for Poly shapes
class PolyShape implements Shape {

    private Color color;
    private Poly poly;

    // Creates a new random polygon
    public PolyShape(int width, int height) {
        color = new Color(Rand.nextFloat(), Rand.nextFloat(),
                Rand.nextFloat(), Rand.nextFloat());
        poly = new Poly(width, height);
    }

    // Creates a new polygon with a specific color and geometry
    public PolyShape(Color color, Poly poly) {
        this.color = color;
        this.poly = poly;
    }

    public PolyShape(PolyShape p) {
        this.color = new Color(p.color.getRed(), p.color.getGreen(), p.color.getBlue(), p.color.getAlpha());
        this.poly = new Poly(p.poly);
    }

    // Draws the polygon onto a Graphics object
    public void draw(Graphics g, int scale) {
        g.setColor(color);
        try {
            g.fillPolygon(poly.getPolygon(scale));
        } catch (Exception e) {
            System.out.println("poly npoints " + poly.getPolygon(scale).npoints + " xpoints " + poly.getPolygon(scale).xpoints.length);
            System.exit(0);
        }
    }

    // Returns a new mutated polygon
    public Shape mutate() {
        // This code returns a new polygon with the same shape and color
        // Change this to mutate the color or shape
        int[] comp = new int[4];

        comp[0] = color.getRed();
        comp[1] = color.getGreen();
        comp[2] = color.getBlue();
        comp[3] = color.getAlpha();

        int redMutationRate = 750;
        int greenMutationRate = 700;
        int blueMutationRate = 550;
        int alphaMutationRate = 500;
        int polyMutationRate = 750;

        if (Rand.nextInt(redMutationRate) == 1) {
            comp[0] = Rand.nextInt(255);
        }
        if (Rand.nextInt(greenMutationRate) == 1) {
            comp[1] = Rand.nextInt(255);
        }
        if (Rand.nextInt(blueMutationRate) == 1) {
            comp[2] = Rand.nextInt(255);
        }
        if (Rand.nextInt(alphaMutationRate) == 1) {
            comp[3] = Rand.nextInt(255);
        }

        if (Rand.nextInt(polyMutationRate) == 1) {
            poly.mutate();
        }


        return new PolyShape(new Color(comp[0], comp[1], comp[2], comp[3]),
                poly);
    }
}

// This is the interface to a factory that creates shapes
// You can add another shape factory if you add a class that represents a
// new type of shape
// Factory class for Poly shapes
interface ShapeFactory {

    public Shape newShape(int width, int height);

    public Shape newShape(Shape get);
}

// Factory class for Poly shapes
// Represents an approximation of the original image using shapes
class PolyFactory implements ShapeFactory {

    public Shape newShape(int width, int height) {
        return new PolyShape(width, height);
    }

    @Override
    public Shape newShape(Shape get) {
        return new PolyShape((PolyShape) get);
    }
}

// Represents an approximation of the original image using shapes
class ImageApprox implements Comparable {

    private int p;
    private final int minPoly = 100;
    private final int maxPoly = 100;
    private BufferedImage original;
    private ShapeFactory factory;
    private int width, height;
    private double cachedFitness;
    // This is the genome -- an ordered list of shapes.
    private List<Shape> shapes;

    // Makes a random image approximation
    public ImageApprox(BufferedImage original, ShapeFactory factory) {
        this.original = original;
        this.factory = factory;
        width = original.getWidth();
        height = original.getHeight();
        cachedFitness = Double.NaN;
        shapes = new ArrayList<Shape>();

        p = minPoly; //10 polygons initially

        for (int i = 0; i < p; i++) {
            shapes.add(factory.newShape(width, height));
        }

    }

    private ImageApprox(ImageApprox aThis) {
        this.original = aThis.original;
        this.factory = aThis.factory;
        this.width = aThis.width;
        this.height = aThis.height;
        this.cachedFitness = aThis.cachedFitness;
        this.shapes = new ArrayList<Shape>();
        p = aThis.p;

        for (int i = 0; i < p; i++) {
            shapes.add(factory.newShape(aThis.shapes.get(i)));
        }

    }

    //Returns the fitness value
    public double getFitness() {
        return cachedFitness;
    }

    // Returns the rendered image approximation
    public BufferedImage image() {
        return image(1);
    }

    // Returns the rendered image approximation enlarged by scale
    public BufferedImage image(int scale) {
        BufferedImage i = new BufferedImage(scale * width, scale * height,
                BufferedImage.TYPE_INT_RGB);
        Graphics g = i.getGraphics();
        for (Shape shape : shapes) {
            shape.draw(g, scale);
        }
        return i;
    }

    // Returns the fitness of the individual
    public double fitness() {
        // This code finds the red, green, and blue components of each pixel of
        // the original image and the approximated image, but always returns a
        // fitness of 0
        // Change this to implement the fitness calculation specified in the
        // assignment
        if (!Double.isNaN(cachedFitness)) {
            return cachedFitness;
        }
        BufferedImage approx = image();
        int[] origPixels = original.getData().
                getPixels(0, 0, original.getWidth(),
                original.getHeight(), (int[]) null);
        int[] apprPixels = approx.getData().
                getPixels(0, 0, original.getWidth(),
                original.getHeight(), (int[]) null);
        long targetDiff = 0, s = 0, t = 0, temp;
        for (int i = 0; i < origPixels.length; i++) {

            temp = origPixels[i] - apprPixels[i];
            temp *= temp;
            targetDiff += temp;


            temp = origPixels[i] * origPixels[i];
            t += temp;

            temp = apprPixels[i] * apprPixels[i];
            s += temp;

        }

        approx = null;

        double d_targetDiff = Math.sqrt(targetDiff);
        double d_s = Math.sqrt(s);
        double d_t = Math.sqrt(t);

        cachedFitness = -Math.log(d_targetDiff / (d_s + d_t));
        return cachedFitness;
    }

    // The compareTo methods are used for sorting the individuals by fitness
    @Override
    public int compareTo(Object that) {
        return compareTo((ImageApprox) that);
    }

    public int compareTo(ImageApprox that) {
        return (int) signum(that.fitness() - this.fitness());
    }

    // Mutates the image approximation
    public void mutate() {
        // This code mutates a random shape
        // You can change this to add other mutations


        //Add Polygon        
        int addPolyRate = 100;
        int randNum = new Random().nextInt(addPolyRate);

        if (p <= maxPoly && randNum == 1) {
            int randIndex = new Random().nextInt(p);
            shapes.add(randIndex, factory.newShape(width, height));
            p++;
        }

        //Remove Polygon
        int removePolyRate = 500;
        randNum = new Random().nextInt(removePolyRate);

        if (p >= minPoly && randNum == 1) {
            int randIndex = new Random().nextInt(p);
            shapes.remove(randIndex);
            p--;
        }

        //Move Polygon
        int movePolygonRate = 250;
        randNum = new Random().nextInt(movePolygonRate);
        if (randNum == 1) {
            int randIndex1 = new Random().nextInt(p);
            Shape get1 = shapes.get(randIndex1);
            int randIndex2 = new Random().nextInt(p);
            Shape get2 = shapes.get(randIndex2);
            shapes.set(randIndex1, get1);
            shapes.set(randIndex2, get2);

        }

        //Mutate color and points of each polygon
        for (int i = 0; i < shapes.size(); i++) {
            int n = Rand.nextInt(shapes.size() - 1);
            shapes.set(n, shapes.get(n).mutate());
        }
        cachedFitness = Double.NaN;
    }

    public List<Shape> GetShapes() {
        return shapes;
    }

    ImageApprox crossover(ImageApprox parent2) throws CloneNotSupportedException {
        ImageApprox child = new ImageApprox(this);

        List<Shape> shapesParent2 = parent2.GetShapes();
        int maxRandCOPoint = p < shapesParent2.size() ? p : shapesParent2.size();
        int randCrossoverPt = Rand.nextInt(maxRandCOPoint);

        child.DoCrossover(randCrossoverPt, shapesParent2);

        return child;
    }

    private void DoCrossover(int randCrossoverPt, List<Shape> shapesParent2) {

        int minSize = p < shapesParent2.size() ? p : shapesParent2.size();

        int maxSize = p > shapesParent2.size() ? p : shapesParent2.size();

        for (int i = randCrossoverPt; i < maxSize; i++) {
            if (p < shapesParent2.size() && i >= p) {
                shapes.add(i, new PolyShape((PolyShape) shapesParent2.get(i)));
            }
            if (p < shapesParent2.size() && i < p) {
                shapes.set(i, new PolyShape((PolyShape) shapesParent2.get(i)));
            }
            if (p > shapesParent2.size() && i < shapesParent2.size()) {
                shapes.set(i, new PolyShape((PolyShape) shapesParent2.get(i)));
            }
        }

        p = p > shapesParent2.size() ? p : shapesParent2.size();

    }
}

class SortedArrayList<T> extends ArrayList<T> {

    @SuppressWarnings("unchecked")
    public void insertSorted(T value) {
        add(value);
        Comparable<T> cmp = (Comparable<T>) value;
        for (int i = size() - 1; i > 0 && cmp.compareTo(get(i - 1)) < 0; i--) {
            T tmp = get(i);
            set(i, get(i - 1));
            set(i - 1, tmp);
        }
    }

    public void removeRange(int start, int end) {
        super.removeRange(start, end);
    }
}

// Displays the original image and the best current approximation
// Displays the original image and the best current approximation
class GeneticSearch {

    private BufferedImage original;
    private int n, k, t;
    private SortedArrayList<ImageApprox> population;
    private int populationSize;

    // Creates a population of n random individuals
    public GeneticSearch(BufferedImage original, int n, int k) {
        this.original = original;
        this.n = n; // population size
        this.k = k; // number of children per generation
        t = 0; //number of generations
        //  population = new ImageApprox[n];
        population = new SortedArrayList<ImageApprox>();
        populationSize = n * 5;
        //n*100 is initial population size
        for (int i = 0; i < populationSize; i++) {
            population.insertSorted(new ImageApprox(original, new PolyFactory()));
        }
    }

    // Makes a new generation of individuals and returns the effort to do so
    // The effort is equal to the number of children in the generation
    public int step() {

        t++;

        //int rand_pop_size = new Random().nextInt(n);

        List<ImageApprox> parent = new ArrayList<ImageApprox>(n);
        List<ImageApprox> children = new ArrayList<ImageApprox>(k);

        double totalFitness = 0.0;

        //Roulette Wheel Selection
        for (ImageApprox i : population) {
            totalFitness += i.getFitness();
        }

        TreeMap cumulProbOfPopulationMap = new TreeMap();
        double cumulProb = 0.0;

        for (ImageApprox i : population) {
            cumulProb += i.getFitness() / totalFitness;
            double currentCumulProb = i.getFitness() / totalFitness;
            cumulProbOfPopulationMap.put(cumulProb - currentCumulProb, i);
        }


        //Choose n parents from population
        for (int i = 0; i < n; i++) {


            //Choose a parent
            double randParent = new Random().nextDouble();

            Map.Entry e = cumulProbOfPopulationMap.floorEntry(randParent);
            ImageApprox parentChosen = (ImageApprox) e.getValue();

            parent.add(parentChosen);

        }
        //Do k crossovers of n parents

        for (int i = 0; i < k; i++) {
            try {
                //Choose two parents randomly
                int parent1Idx = new Random().nextInt(n);
                int parent2Idx = new Random().nextInt(n);

                ImageApprox parent1 = parent.get(parent1Idx);
                ImageApprox parent2 = parent.get(parent2Idx);

                ImageApprox child = parent1.crossover(parent2);

                child.mutate();

                children.add(child);

            } catch (CloneNotSupportedException ex) {
                Logger.getLogger(GeneticSearch.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        //Insert children to population

        for (int i = 0; i < children.size(); i++) {
            population.insertSorted(children.get(i));
        }

        //Now trim population to original size

        int currentPopSize = population.size();
        population.removeRange(populationSize, currentPopSize);

        return n + t * k;
    }

    // Returns the original image
    public BufferedImage getOriginal() {
        return original;
    }

    // Returns the best current image approximation, assuming the individuals
    // are sorted in decreasing order of fitness.
    public ImageApprox getBest() {
        return population.get(0);
    }
}

// Displays the original image and the best current approximation
class Visualization extends Frame {

    public static final long serialVersionUID = 2011091201L;
    private static int TopBorder = 32;
    private GeneticSearch s;
    private int scale;

    public Visualization(GeneticSearch s) {
        super("Visualization");
        this.s = s;
        scale = 1;
        int w = s.getOriginal().getWidth();
        int h = s.getOriginal().getHeight();
        while ((scale + 1) * w <= 600 && (scale + 1) * h <= 800) {
            scale++;
        }
        setSize(2 * scale * w, scale * h + TopBorder);
        setVisible(true);
    }

    public void paint(Graphics g) {
        Dimension dim = getSize();
        int w = s.getOriginal().getWidth();
        int h = s.getOriginal().getHeight();
        g.drawImage(s.getOriginal().
                getScaledInstance(scale * w, scale * h, Image.SCALE_FAST),
                0, TopBorder, this);
        g.drawImage(s.getBest().image(scale), scale * w, TopBorder, this);
    }

    public void update(Graphics g) {
        paint(g);
    }
}

public class A1Main {

    public static void main(String[] args) throws IOException {
        final int reportFreq = 1000;
        int n = 0, k = 0;
        long e = 0;
        String filename = null;
        try {
            n = Integer.parseInt(args[0]);
            k = Integer.parseInt(args[1]);
            e = Long.parseLong(args[2]);
            filename = args[3];
        } catch (Exception ex) {
            System.out.println("Usage: java A1Main <n> <k> <e> <image>");
            System.exit(0);
        }
        try {


            BufferedImage image = ImageIO.read(new File(filename));
            GeneticSearch s = new GeneticSearch(image, n, k);
            Visualization v = new Visualization(s);

            //Add label

            Label l = new Label("Generation");
            Label l1 = new Label("                    ");


            v.setLayout(new FlowLayout(FlowLayout.CENTER));

            v.add(l);
            v.add(l1);

            Apng apng = ApngFactory.createApng();
            apng.setPlayCount(1);
            long effort = 0;
            int effortReport = reportFreq;
            int o = 0;
            while (effort <= e) {
                v.repaint();
                if (effort >= effortReport) {
                    apng.addFrame(s.getBest().image(), 1000 / 15);
                    effortReport += reportFreq;
                }

                effort += s.step();

                l1.setText(Long.toString(effort));
            }

            System.out.println("After " + effort + " generations, Highest fitness value is " + s.getBest().fitness());
            apng.assemble(new File(filename + ".png"));
            v.setVisible(false);



        } catch (IOException ioe) {
            System.out.println("IO Error");
        } catch (JapngException je) {
            System.out.println("Japng Error");
        }

        System.exit(0);
    }
}
