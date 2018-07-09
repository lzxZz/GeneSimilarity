import java.util.ArrayList;



//保存注释项的数据结构
public class AnnoNode {
    public String               GeneName;
    public String               NameSpace;
    public String               EvidenceCode;
    public String               Type;   //绝大多数都是gene
    public String               ID;
    public ArrayList<String>    Synonym;    //保存同义词，同义词为ID
    public ArrayList<String>    RefIDS;


    public boolean isSynonym(String name){
        if (name.equals(this.GeneName))
            return true;

        for (String n:Synonym){
            if (name.equals(n))
                return  true;
        }

        return false;
    }

    public AnnoNode(){
       // OntoIds = new ArrayList<String>();
        Synonym = new ArrayList<String>();
        RefIDS= new ArrayList<String>();
    }
    public AnnoNode(String id){
       // OntoIds = new ArrayList<String>();
        Synonym = new ArrayList<String>();
        RefIDS= new ArrayList<String>();
    }
    public AnnoNode(String name,String id){
      //  OntoIds = new ArrayList<String>();
        Synonym = new ArrayList<String>();
      //  OntoIds.add((id));
        this.ID  = id;
        this.GeneName = name;
        RefIDS= new ArrayList<String>();
    }
}
