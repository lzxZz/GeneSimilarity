import java.lang.reflect.Array;
import java.util.ArrayList;

//术语项数据结构，保存了术语的名称，ID，所属分支，定义，是否过时，以及父节点ID列表和父节点节点列表
public class TermNode {
    public String                   ID
;   public String                   name
;   public String                   namespace
;   public String                   def
;   public boolean                  isObsolete
;   public ArrayList<TermNode>      IParentNodes
;   public ArrayList<TermNode>      PParentNodes
;   public ArrayList<String>        IParentIDs
;   public ArrayList<String>        PParentIDs
;

    public TermNode(){
        isObsolete      = false;
;       IParentIDs      = new ArrayList<String>()
;       PParentIDs      = new ArrayList<String>()
;       IParentNodes    = new ArrayList<TermNode>()
;       PParentNodes    = new ArrayList<TermNode>()
;
    }
}
