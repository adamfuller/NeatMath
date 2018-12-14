
/**
 *  0.0.5   12/13/2018  Added Curve class with value(double) value(int) derivative and toString methods
 *  0.0.4   12/13/2018  Added NoSolutionException
 *  0.0.3   12/13/2018  Added pointOfIntersection2D methods based on 4 points and 2 points and 2 vectors
 *  0.0.2   12/13/2018  Added Point and Vector Class
 *  0.0.1   3/12/2018   Original Coding
 */
import java.awt.Color;

/**
 * Class containing cool math tools I have needed
 */
public class AdamMath {
    /*
     * Returns the distance between two x,y coordinates
     */
    public double distance(double x1, double y1, double x2, double y2) {
        return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    }

    /*
     * Returns the square of the distance between two x,y coordinates
     */
    public double distanceSq(double x1, double y1, double x2, double y2) {
        return (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
    }

    /*
     * Returns the average value between two (x,y) returns format {xValue, yValue}
     */
    public double[] midPoint2D(double x1, double y1, double x2, double y2) {
        double[] returnVal = new double[2];
        returnVal[0] = (x1 + x2) / 2;
        returnVal[1] = (y1 + y2) / 2;

        return returnVal;
    }

    /*
     * Returns the average value between two arrays of (x,y) coordinates returns
     * format {xMidPoint, yMidPoint}
     */
    public double[] midPoint2D(double[] x, double[] y) {
        double[] returnVal = new double[2];
        if (x.length == y.length) {
            for (int i = 0; i < x.length; i++) {
                returnVal[0] += x[i];
                returnVal[1] += y[i];
            }
            returnVal[0] /= x.length;
            returnVal[1] /= y.length;
        }
        return returnVal;
    }

    /*
     * finds and returns the closest point between a point and a line segment in the
     * format {pointX, pointY}
     */
    public double[] closestPoint2D(double lineX1, double lineY1, double lineX2, double lineY2, double pointX,
            double pointY) {
        double[] returnValue = new double[2];
        // calculates
        double A1 = lineY2 - lineY1;
        double B1 = lineX1 - lineX2;
        double C1 = (lineY2 - lineY1) * lineX1 + (lineX1 - lineX2) * lineY1;
        double C2 = -B1 * pointX + A1 * pointY;
        double det = A1 * A1 - -B1 * B1;

        if (det != 0) {
            returnValue[0] = ((A1 * C1 - B1 * C2) / det);
            returnValue[1] = ((A1 * C2 - -B1 * C1) / det);
        } else {
            returnValue[0] = pointX;
            returnValue[1] = pointY;
        }
        return returnValue;
    }

    /**
     * Returns a point of intersection from vectors from points a to b and c to d.
     * <p>
     * Points a and b are reversable without affecting the outcome as well as c and
     * d. Returns null if the points are parallel
     * 
     * @param a first point
     * @param b second point
     * @param c third point
     * @param d fourth point
     * @return {@code Point} object of the intersection point
     * @throws NoSolutionException if the points produce parallel vectors
     */
    public Point pointOfIntersect2D(Point a, Point b, Point c, Point d) throws NoSolutionException {
        Vector ab = new Vector(a, b);
        Vector cd = new Vector(c, d);

        if (ab.equals(cd)) { // vectors are the same so a->b and c->d are parallel
            throw new NoSolutionException("Input points produce parallel vectors");
        }

        double xInt = ((c.y - c.x * (cd.y / cd.x)) - (b.y - b.x * (ab.y / ab.x))) / ((ab.y / ab.x) - (cd.y / cd.x));
        double yInt = xInt * (ab.y / ab.x) + (b.y - b.x * (ab.y / ab.x));

        return new Point(xInt, yInt);
    }

    /**
     * Returns the point of intersection based on vectors aVec and bVec that
     * correlate with points a and b respectively.
     * <p>
     * 
     * @param a    point aVec originates from
     * @param b    point bVec originates from
     * @param aVec direction of ray from point a
     * @param bVec direction of ray from point b
     * @return {@code Point} object of intersection or {@code null} if parallel
     * @throws NoSolutionException if the Vectors are parallel or rays don't
     *                             intersect
     */
    public Point pointOfIntersect2D(Point a, Point b, Vector aVec, Vector bVec) throws NoSolutionException {
        if (aVec.equals(bVec)) { // vectors are the same so rays are parallel
            throw new NoSolutionException("Input vectors are parallel");
        }

        // Following only valid assuming they intersect
        double xInt = ((a.y - a.x * (aVec.y / aVec.x)) - (b.y - b.x * (bVec.y / bVec.x)))
                / ((bVec.y / bVec.x) - (aVec.y / aVec.x));
        double yInt = xInt * (bVec.y / bVec.x) + (b.y - b.x * (bVec.y / bVec.x));

        Point intPoint = new Point(xInt, yInt);

        // unit vector of a->intPoint should equal unit vector of aVec
        // and unit vector of b->intPoint should equal unit vector of bVec
        // if the rays never intersect throw an exception
        if (!(aVec.unit().equals(new Vector(a, intPoint).unit())
                && bVec.unit().equals(new Vector(b, intPoint).unit()))) {
            throw new NoSolutionException("Vectors never intersect");
        }

        return intPoint;
    }

    /*
     * Cosine approximation using Mclaurin series
     */
    private double cos(double angle) {
        double ans = 0.0;
        int sign = 1;
        for (int i = 0; i < 20; i += 1) {
            ans += sign * (double) (power(angle, 2 * i) / factorial(2 * i));
            sign *= -1;
        }

        return ans;
    }

    /**
     * Returns a linear interpolation between two numbers, {@code minNum} and
     * {@code maxNum}, to two other numbers {@code minMap} and {@code maxMap}
     * 
     * @param num
     * @param minNum
     * @param maxNum
     * @param minMap
     * @param maxMap
     * @return
     */
    public static double map(double num, double minNum, double maxNum, double minMap, double maxMap) {
        /*
         * double slope = (maxMap - minMap)/(maxNum - minNum); double dif1 = (num -
         * minNum); double move = slope*dif1;
         */
        double allInOne = ((maxMap - minMap) / (maxNum - minNum)) * (num - minNum) + minMap;
        // y = (y2 - y1)/(x2 - x1) * (x-x1) + y1
        return (allInOne);
    }

    /**
     * Returns color based on closeness to min and max
     * 
     * @param val
     * @param min
     * @param max
     * @return
     */
    public Color getColor(int val, int min, int max) {
        int r = 0;
        int g = 0;
        int b = 0;
        int mid = (max + min) / 2;

        if (val > mid) {
            r = (int) map(val * 1.0, mid * 1.0, max * 1.0, 0, 254);
            g = (int) map(val * 1.0, mid * 1.0, max * 1.0, 254, 0);
        } else if (val < mid) {
            g = (int) map(val * 1.0, min * 1.0, mid * 1.0, 0, 254);
            b = (int) map(val * 1.0, min * 1.0, mid * 1.0, 254, 0);
        } else {
            g = 254;
        }
        return new Color(r, g, b);
    }

    /*
     * Sine approximation using Mclaurin Series
     */
    private double sin(double angle) {
        double ans = 0.0;
        int sign = 1;
        for (int i = 0; i < 10; i++) {
            ans += sign * (double) (power(angle, 2 * i + 1) / factorial(2 * i + 1));
            sign *= -1;
        }

        return ans;
    }

    /*
     * Returns num^exp
     */
    public static double power(double num, int exp) {
        double ans = 1;
        for (int i = 0; i < exp; i++) {
            ans *= num;
        }
        return ans;
    }

    /*
     * Returns num! ie. 3! = 3*2*1
     */
    public long factorial(int num) {
        long ans = 1;
        for (int i = num; i > 0; i--) {
            ans *= i;
        }
        return ans;
    }

    /*
     * Returns a weighted average of a 1D array
     */
    public static double[] weighted(double[] array) {
        // create duplicate of the same size
        double[] weighted = new double[array.length];
        // sum of the whole input array
        double sum = 0;
        // weight of the current position
        double weight = 1.0;
        // change in the weight
        double deltaWeight = 1.0 / (array.length - 1);
        // total amount of weight for each position
        double totalWeighting = 0.0;

        for (int i = 0; i < array.length; i++) {
            sum += array[i];
        }

        // iterate through each element in the weighted array
        for (int i = 0; i < weighted.length; i++) {
            totalWeighting = 0.0;

            for (int n = 0; n < weighted.length; n++) {
                // sum up weighting for this element
                totalWeighting += 1.0 - (1.0 / (array.length - 1) * abs(i - n));

            }
            // System.out.println(totalWeighting);

            // iterate through each element in the input array
            for (int j = 0; j < array.length; j++) {

                // calculate weight at the position
                weight = 1.0 - (deltaWeight * abs((i - j)));

                System.out.println(weight / totalWeighting);
                weighted[i] += (array[j] * weight) / totalWeighting;
            }
            // weighted[i] /= sum / (array.length);
        }

        return weighted;

    }

    /*
     * Returns the absolute value of a double input
     */
    public static double abs(double num) {
        if (num >= 0) {
            return num;
        }
        return num * -1.0;
    }

    /*
     * Returns the absolute value of an int input
     */
    public static int abs(int num) {
        if (num >= 0) {
            return num;
        }
        return num * -1;
    }

    /*
     * Returns square root approximation of a double input
     */
    public static double sqrt(double num) {
        double leftPoint = 0.0;
        double midPoint = num / 2.0;
        double rightPoint = num;
        double leftMidPoint = (leftPoint + midPoint) / 2.0;
        double rightMidPoint = (rightPoint + midPoint) / 2.0;
        int iter = 0;

        if (num <= 0.000001) {
            return 0.0;
        }

        // repeat untill desired accuracy or iteration max
        while (abs(num - midPoint * midPoint) >= 0.0000000001 && iter < 100) {
            leftMidPoint = (leftPoint + midPoint) / 2.0;
            rightMidPoint = (rightPoint + midPoint) / 2.0;
            // left mid point is closer to square root
            if (abs(num - leftMidPoint * leftMidPoint) < abs(num - rightMidPoint * rightMidPoint)) {
                rightPoint = midPoint;
                midPoint = leftMidPoint;

            } else {
                // right mid point is closer to square root
                leftPoint = midPoint;

                // switch the midpoint to the midpoint between the current midpoint and the
                // right point
                midPoint = rightMidPoint;
            }

            // boost iteration counter
            iter++;

            // System.out.printf("%d,%f,%f,%f\n",iter, leftPoint, midPoint, rightPoint);
        }
        return midPoint;
    }

    /*
     * Returns the natural log of an input number
     */
    public static double ln(double num) {
        num -= 1.0;
        double approx = 0.0;
        int sign = 1;
        for (int i = 1; i < 25; i++) {
            approx += sign * power(num, i) / i;
            sign *= -1;
        }
        return approx;
    }

    public static void main(String[] args) {
        Curve c = new Curve(new double[] { 1, -2, 1 });

        System.out.println(c.unit(1.0).toString());
        System.out.println(c.toString());

    }

    public static class Point {
        public double x, y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        public Point(int x, int y) {
            this.x = (double) x;
            this.y = (double) y;
        }

    }

    public static class Vector {
        public double x, y;

        @Override
        public boolean equals(Object otherVector) {
            if (otherVector.getClass() != Vector.class) {
                return false;
            }
            Vector otherVector_ = (Vector) otherVector;
            if (Math.abs(this.x - otherVector_.x) < 0.000000001) {
                if (Math.abs(this.y - otherVector_.y) < 0.000000001) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public String toString(){
            String result = "";
            result+=String.valueOf(this.x) + "x";
            result+= (this.y>0?"+":"")+String.valueOf(this.y)+"y";
            return result;
        }

        /**
         * Create new Vector object from integer x and y values
         * 
         * @param x
         * @param y
         */
        public Vector(int x, int y) {
            this.x = (double) x;
            this.y = (double) y;
        }

        /**
         * Create new Vector object from double x and y values
         * 
         * @param x
         * @param y
         */
        public Vector(double x, double y) {
            this.x = x;
            this.y = y;
        }

        /**
         * Creates a vector object from Point {@code a} to Point {@code b}
         * 
         * @param a initial point
         * @param b final point
         */
        public Vector(Point a, Point b) {
            this.x = b.x - a.x;
            this.y = b.y - a.y;
        }

        public Vector unit() {
            double mag = Math.sqrt(this.x * this.x + this.y * this.y);
            return new Vector((this.x / mag), (this.y / mag));
        }
    }

    public static class Curve {
        private double coeff[]; // coefficients of the curve {1,2,3} = 3x^2+2x+1

        public Curve(double coeff[]) {
            this.coeff = coeff;
        }

        public Curve(String equation){
            
        }

        @Override
        public String toString(){
            StringBuilder equation = new StringBuilder();

            for (int power = this.coeff.length-1; power>=0; power--){
                if (this.coeff[power] > 0){
                    if (equation.length() != 0){
                        equation.append("+");
                    }
                    equation.append(this.coeff[power]);
                    if (power >= 1 ){
                        equation.append("x");
                        if (power > 1){
                            equation.append("^");
                            equation.append(power);
                        }
                    }
                } else if (this.coeff[power] < 0){
                    equation.append(this.coeff[power]);
                    if (power >= 1){
                        equation.append("x");
                        if (power > 1){
                            equation.append("^");
                            equation.append(power);
                        }
                    }

                }
            }

            return equation.toString();
        }

        /**
         * Gets the value of the curve at {@code x}
         * @param x - independent variable of the curve
         * @return {@code double} value at the input value
         */
        public double value(double x) {
            double val = 0.0;
            for (int power = 0; power < this.coeff.length; power++) {
                val += this.coeff[power] * AdamMath.power(x, power);
            }
            return val;
        }

        /**
         * Gets the value of the curve at {@code x}
         * @param x - independent variable of the curve
         * @return {@code double} value at the input value
         */
        public double value(int x) {
            double val = 0.0;
            for (int power = 0; power < this.coeff.length; power++) {
                val += this.coeff[power] * AdamMath.power(x, power);
            }
            return val;
        }

        /**
         * Creates a curve of the derivative of this curve
         * @return Curve that is the derivative of this curve
         */
        public Curve derivative(){
            if (this.coeff.length == 1){
                return new Curve(new double[]{0.0});
            }
            double newCoeff[] = new double[this.coeff.length-1];

            for (int power = 1; power < coeff.length; power++){
                newCoeff[power-1] = this.coeff[power]*power;
            }

            return new Curve(newCoeff);
        }

        /**
         * Returns the tangent unit vector at x on this curve
         * <p>
         * The vector is pointing from left to right (-x to +x) if applicable
         * @param x - position of the tangent vector
         * @return {@code Vector} object tangent to this curve at {@code x}
         */
        public Vector unit(double x){
            Curve dCurve = this.derivative();
            double slope = dCurve.value(x);

            double mag = Math.sqrt(1+slope*slope);

            return new Vector(1/mag, slope/mag);
        }
    }

    /**
     * Exception when inputs produce no solution for the chosen function
     */
    private class NoSolutionException extends Exception {
        public NoSolutionException(String message) {
            super(message);
        }
    }

}