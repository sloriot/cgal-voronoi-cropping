import CGAL.Kernel.Point_2;
import CGAL.Kernel.Iso_rectangle_2;
import CGAL.Voronoi_cropping_2.HalfedgeDS;
import CGAL.Voronoi_cropping_2.HDS_Face_handle;
import CGAL.Voronoi_cropping_2.HDS_Halfedge_handle;
import CGAL.Voronoi_cropping_2.CGAL_Voronoi_cropping_2;
import CGAL.Voronoi_cropping_2.Voronoi_cropping_2;
import java.util.Vector;

public class test_vcrop {

  public static void print_info(HalfedgeDS hds)
  {
    System.out.println("HDS-info "+hds.size_of_vertices()+", "+hds.size_of_halfedges()+", "+hds.size_of_faces());
  }


  public static void main(String arg[]){

    Vector<Point_2> points = new Vector<Point_2>(7);
    points.add( new Point_2(0,0) );
    points.add( new Point_2(0,1) );
    points.add( new Point_2(1,1) );

    int[] colors = new int[] {1,1,2};

    Iso_rectangle_2 isorect =
      new Iso_rectangle_2( new Point_2(-2,-2), new Point_2(2,2) );

    HalfedgeDS hds = new HalfedgeDS();


    Voronoi_cropping_2 vcrop=new Voronoi_cropping_2();
    vcrop.insert(points.iterator(), colors);

    vcrop.voronoi_diagram(hds, isorect);
    print_info(hds);
    CGAL_Voronoi_cropping_2.join_faces(hds);
    print_info(hds);

    vcrop.voronoi_diagram(hds);
    print_info(hds);
    CGAL_Voronoi_cropping_2.join_faces(hds);
    print_info(hds);

    for (HDS_Face_handle fh : hds.faces() )
    {
      if (fh.color() == 2) fh.set_color(1); //update the color
    }

    vcrop.insert( new Point_2(1,0), 1);
    vcrop.insert( new Point_2(0.5,0.5), 2);
    vcrop.voronoi_diagram(hds, isorect);//clear the hds and replace it with a new
    CGAL_Voronoi_cropping_2.join_faces(hds);
    print_info(hds);



    for ( HDS_Face_handle fh : hds.faces() )
    {
      if ( fh.has_holes() ) System.out.println("Start a new face with holes");
      else System.out.println("Start a new face");
      System.out.println("Corresponding input point: "+vcrop.get_point(fh.point_index())+"\n");
      HDS_Halfedge_handle hedge = fh.halfedge();
      HDS_Halfedge_handle first = hedge.clone();

      do{
          System.out.println("2 "+ hedge.vertex().point() + " 0 " + hedge.opposite().vertex().point() + " 0");
          hedge=hedge.next();
      }
      while(!hedge.equals(first) );

      if ( fh.has_holes() )
      {
        for (HDS_Halfedge_handle hh : fh.holes() )
        {
          first = hh.clone();

          do{
              System.out.println("2 "+ hh.vertex().point() + " 0 " + hh.opposite().vertex().point() + " 0");
              hh=hh.next();
          }
          while(!hh.equals(first) );
        }
      }
    }

  }

}
