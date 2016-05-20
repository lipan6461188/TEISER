typedef struct _ENTRY
{
  unsigned int used;
  ENTRY entry;  //定义了key和data的结构体
} _ENTRY ;


struct my_hsearch_data
{
  int filled;
  _ENTRY* table;  //放置任意多个key-value
  int size;
}; 

int  my_hcreate_r  (size_t nel,  struct my_hsearch_data *htab);
//仅仅是把哈希表中的table释放掉，而不释放哈希对象本身
void my_hdestroy_r (struct my_hsearch_data *htab);
int  my_hsearch_r  (ENTRY item, ACTION action, ENTRY **retval, struct my_hsearch_data *htab);
     //ACTION中定义了是输入还是查找

