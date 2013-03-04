import CGAL.Kernel.Point_2;
import CGAL.Kernel.Iso_rectangle_2;
import CGAL.Voronoi_cropping_2.HalfedgeDS;
import CGAL.Voronoi_cropping_2.HDS_Face_handle;
import CGAL.Voronoi_cropping_2.HDS_Halfedge_handle;
import CGAL.Voronoi_cropping_2.CGAL_Voronoi_cropping_2;
import java.util.Vector;

public class test_vcrop {

  public static void main(String arg[]){
    
    Vector<Point_2> points = new Vector<Point_2>(7);
    points.add( new Point_2(0,0) );
    points.add( new Point_2(0,1) );
    points.add( new Point_2(1,1) );

    int[] colors = new int[] {1,1,2};
    
    Iso_rectangle_2 isorect =
      new Iso_rectangle_2( new Point_2(-2,-2), new Point_2(2,2) );
    
    HalfedgeDS hds = new HalfedgeDS();
    
    CGAL_Voronoi_cropping_2.cropped_voronoi_diagram_2(points.iterator(), colors, isorect, hds);
    CGAL_Voronoi_cropping_2.join_faces(hds);
   
    for ( HDS_Face_handle fh : hds.faces() )
    {
      HDS_Halfedge_handle hedge = fh.halfedge();
      HDS_Halfedge_handle first = hedge.clone();
      
      do{
          System.out.println("2 "+ hedge.vertex().point() + " 0 " + hedge.opposite().vertex().point() + " 0");
          hedge=hedge.next();
      }
      while(!hedge.equals(first) );
    }
    
  }
  
}
