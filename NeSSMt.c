#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<time.h>
#include<sys/stat.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define N 500 //max replicate
int FileSize(char *file)
{
	struct stat statbuf;
	stat(file,&statbuf);
	return statbuf.st_size;
}

void minus_strand(char *buffle,char *minus_buffle,int len)
{
	int circle_f;
	for(circle_f=0;circle_f<len;circle_f++)
	{
		if(buffle[circle_f]=='A')
			minus_buffle[len-1-circle_f]='T';
		else if (buffle[circle_f]=='T')
			minus_buffle[len-1-circle_f]='A';
		else if (buffle[circle_f]=='C')
			minus_buffle[len-1-circle_f]='G';
		else
			minus_buffle[len-1-circle_f]='C';
	}
}

struct StoreGene
{
	char *Genome;
	char *Gene;
	char *Fasta;
	int *countTreat;
	int *countControl;
	float FC;
	struct StoreGene *Next;
};

struct InputGene
{
	char Name[100];
	int count;
	float FC;
	struct InputGene *Next;
};

struct Genome
{
	char Name[400]; //max length of name is 400
	int ReadNum[2];
	float DE[5];
	int mRNA; //total mRNA number
	int rRNA;
	int MaxLength; //max length
	char *genome_path;
	char *transcript_path;
	char *annotation_path;
	int AnnotationType; //1:genbank, 2:gff3
	char *rRNA_path;
	struct InputGene *gene_list;
	struct Genome *Next;
};

int split(char **arr,char *str,char *del)
{
	int i=0;
	char *s=NULL;

	s=strtok(str,del);
	while(s!=NULL)
	{
		i++;
		*arr++=s;
		s=strtok(NULL,del);
	}
	return i;
}

struct Genome *Abundance(char *file,float express_per,float up_expression_per,float down_expression_per,float maxFC,float minFC,int readNum,int Flag_control_treat)
{
	int size,i,j,flag,num,left[2];
	char readIn[500],*part[5];
	FILE *fp;
	float par[5],temp,total[2];
	struct Genome *head,*new,*store;
	
	flag=0;
	total[0]=0;
	par[0]=express_per;
	par[1]=up_expression_per;
	par[2]=down_expression_per;
	par[3]=maxFC;
	par[4]=minFC;
	left[0]=readNum;
	left[1]=readNum;
	size=sizeof(struct Genome);
	memset(readIn,'\0',500*sizeof(char));
	fp=fopen(file,"r");
	if(fp==NULL)
	{
		printf("Error! %s\n",file);
		exit(1);
	}
	i=0;
	flag=0;
	while(fgets(readIn,500,fp)!=NULL)
	{
		if(readIn[0]=='[')
		{
			memset(readIn,'\0',500*sizeof(char));
			continue;
		}
		j=strlen(readIn);
		if(readIn[j-1]=='\n')
			readIn[j-1]='\0';
		num=split(part,readIn,"\t");
		if(Flag_control_treat==1&&num==3) //user input abundance change
			flag=1;
		total[0]=total[0]+atof(part[1]);
		memset(readIn,'\0',500*sizeof(char));
	}
	rewind(fp); //begin again
	if(flag==1) //user input abundance change
	{
		total[1]=0;
		while(fgets(readIn,500,fp)!=NULL)
		{
			if(readIn[0]=='[')
			{
				memset(readIn,'\0',500*sizeof(char));
				continue;
			}
			j=strlen(readIn);
			if(readIn[j-1]=='\n')
				readIn[j-1]='\0';
			num=split(part,readIn,"\t");
			if(num==3)//user input abundance change
				total[1]=total[1]+atof(part[2]);
			else //use before abundance
				total[1]=total[1]+atof(part[1]);
			memset(readIn,'\0',500*sizeof(char));
		}
		rewind(fp);
	}
	while(fgets(readIn,500,fp)!=NULL)
	{
		if(readIn[0]=='[')
		{
			split(part,readIn,"=");
			j=strlen(part[1]);
			part[1][j-1]='\0';
			temp=atof(part[1]);
			if(strcmp(part[0],"[Num_Expression")==0)
				par[0]=temp;
			else if(strcmp(part[0],"[Num_Up_DE")==0)
				par[1]=temp;
			else if(strcmp(part[0],"[Num_Down_DE")==0)
				par[2]=temp;
			else if(strcmp(part[0],"[Max_FC")==0)
				par[3]=temp;
			else if(strcmp(part[0],"[Min_FC")==0)
				par[4]=temp;
			else
				printf("Warning: Can't recognize %s!\n",part[0]);
		}
		else
		{
			j=strlen(readIn);
			if(readIn[j-1]=='\n')
				readIn[j-1]='\0';
			num=split(part,readIn,"\t");
			new=(struct Genome *)malloc(size);
			strcpy(new->Name,part[0]);
			temp=atof(part[1]);
			new->ReadNum[0]=(int)(temp*readNum/total[0]+0.5);
			left[0]=left[0]-new->ReadNum[0];
			if(flag==1&&num==2)
				new->ReadNum[1]=(int)(temp*readNum/total[1]+0.5);
			if(flag==1&&num==3)
			{
				temp=atof(part[2]);
				new->ReadNum[1]=(int)(temp*readNum/total[1]+0.5);
			}
			if(flag==1)
				left[1]=left[1]-new->ReadNum[1];
			new->MaxLength=0; //as flag
			new->genome_path=NULL;
			new->annotation_path=NULL;
			new->rRNA_path=NULL;
			new->transcript_path=NULL;
			new->gene_list=NULL;
			new->Next=NULL;
			for(j=0;j<5;j++)
				new->DE[j]=par[j];
			if(i==0)
			{
				head=new;
				store=new;
				i++;
			}
			else
			{
				store->Next=new;
				store=new;
			}
		}
		memset(readIn,'\0',500*sizeof(char));
	}
	fclose(fp);
	if(Flag_control_treat==1&&flag==0)
		new->ReadNum[1]=-1;  //as a flag, the last change into the first after index
	if(left[0]>=0)
		head->ReadNum[0]=head->ReadNum[0]+left[0];
	else
	{
		new=head;
		while(new)
		{
			if(left[0]==0)
				break;
			if(new->ReadNum[0]+left[0]>=0)
			{
				new->ReadNum[0]=new->ReadNum[0]+left[0];
				break;
			}
			left[0]=left[0]+new->ReadNum[0];
			new->ReadNum[0]=0;
			new=new->Next;
		}
	}
	if(Flag_control_treat==1&&flag==1)
	{
		if(left[1]>=0)
			head->ReadNum[1]=head->ReadNum[1]+left[1];
		else
		{
			new=head;
			while(new)
			{
				if(left[1]==0)
					break;
				if(new->ReadNum[1]+left[1]>=0)
				{
					new->ReadNum[1]=new->ReadNum[1]+left[1];
					break;
				}
				left[1]=left[1]+new->ReadNum[1];
				new->ReadNum[1]=0;
				new=new->Next;
			}
		}
	}
	return head;
}

struct Genome *readIndex(char *file,struct Genome *head)
{
	char *part[10],readIn[3000],*name[2];
	FILE *fp;
	int flag,i,j,num,length;
	struct Genome *node,*newhead,*before;
	
	newhead=NULL;
	fp=fopen(file,"r");
	if(fp==NULL)
	{
		printf("Error! Can't open %s file!\n",file);
		exit(1);
	}
	memset(readIn,'\0',3000*sizeof(char));
	while(fgets(readIn,3000,fp)!=NULL)
	{
		length=strlen(readIn);
		if(readIn[length-1]=='\n')
			readIn[length-1]='\0';
		num=split(part,readIn,"\t");
		i=0;
		j=split(name,part[0],":");
		node=head;
		flag=0;
		while(node)
		{
			if(strcmp(node->Name,name[0])==0)
			{
				flag++;
				break;
			}
			if((j==2)&&(strcmp(node->Name,name[1])==0))
			{
				flag++;
				break;
			}
			i++;
			before=node;
			node=node->Next;
		}
		if(flag==0)
		{
			memset(readIn,'\0',3000*sizeof(char));
			continue;  //not in abundance table
		}
	
		node->MaxLength=atoi(part[1]);
		node->mRNA=atoi(part[2]);
		node->rRNA=atoi(part[3]);
		for(j=4;j<num;j++)
		{
			length=strlen(part[j]);
			if(part[j][0]=='T')
			{
				node->transcript_path=(char *)malloc(length*sizeof(char));
				memset(node->transcript_path,'\0',length*sizeof(char));
				strcpy(node->transcript_path,&part[j][2]);
			}
			else if(part[j][0]=='R')
			{
				node->rRNA_path=(char *)malloc(length*sizeof(char));
				memset(node->rRNA_path,'\0',length*sizeof(char));
				strcpy(node->rRNA_path,&part[j][2]);
			}
			else if(part[j][1]==':')
			{
				node->genome_path=(char *)malloc(length*sizeof(char));
				memset(node->genome_path,'\0',length*sizeof(char));
				strcpy(node->genome_path,&part[j][2]);
			}
			else if(part[j][1]=='B') //genbank
			{
				node->annotation_path=(char *)malloc(length*sizeof(char));
				memset(node->annotation_path,'\0',length*sizeof(char));
				strcpy(node->annotation_path,&part[j][3]);
				node->AnnotationType=1;
			}
			else
			{
				node->annotation_path=(char *)malloc(length*sizeof(char));
				memset(node->annotation_path,'\0',length*sizeof(char));
				strcpy(node->annotation_path,&part[j][3]);
				node->AnnotationType=2;
			}
		}
		
		if(i==0) //this node from head
			head=head->Next;
		else
			before->Next=node->Next;
		node->Next=newhead;
		newhead=node;
		memset(readIn,'0',3000*sizeof(char));
		if(head==NULL)
			break; //over
	}
	fclose(fp);
	if(newhead==NULL)
		return NULL;
	if(head!=NULL)
	{
		node=newhead;
		while(node->Next)
			node=node->Next;
		node->Next=head;
	}
	return newhead;
}
//if the input parameters are not right, use the default parameters
void ExpressNum(struct Genome *head,float express_per,float up_expression_per,float down_expression_per)
{
	struct Genome *node;

	node=head;
	while(node)
	{
		if(node->MaxLength==0)
		{
			node=node->Next;
			continue;
		}
		if(node->DE[0]<1)
			node->DE[0]=node->DE[0]*node->mRNA+0.5;
		if(node->DE[0]<1)
			node->DE[0]=node->mRNA;
		if(node->DE[0]>node->mRNA)
			node->DE[0]=express_per*node->mRNA+0.5;
		if(node->DE[1]<1) //up
			node->DE[1]=node->DE[1]*node->DE[0]+0.5;
		if(node->DE[2]<1)
			node->DE[2]=node->DE[2]*node->DE[0]+0.5;
		if(node->DE[1]+node->DE[2]>node->DE[0]) //diff
		{
			node->DE[1]=node->DE[0]*up_expression_per+0.5;
			node->DE[2]=node->DE[0]*down_expression_per+0.5;
		}
		node=node->Next;
	}
}

int checkFlag(char readIn[],char search[],int flag) //flag=1: check \d+..; flag=0, don't check
{
	int i,size,length;

	size=strlen(search);

	i=0;
	while(readIn[i]==' ')
		i++;
	if(memcmp(&readIn[i],search,size)!=0)
		return 0;

	if(flag==0)
		return 1;

	length=strlen(readIn);
	while(i<length-1)
	{
		if(readIn[i]!='.')
		{
			i++;
			continue;
		}
		if(readIn[i+1]!='.')
		{
			i++;
			continue;
		}
		if(readIn[i-1]>='0'&&readIn[i-1]<='9')
			return 1;
		i++;
	}
	return 0;
}

int getint(char pos[])
{
	int i,length;

	length=strlen(pos);
	i=0;
	while(pos[i]<'0'||pos[i]>'9')
		i++;
	length--;
	while(pos[length]>'9'||pos[length]<'0')
	{
		pos[length]='\0';
		length--;
	}
	return atoi(&pos[i]);
}

int getLength(char PosString[])
{
	int length,start,end,i,num;
	char *part[200],*pos[2],temp[2000];

	memset(temp,'\0',2000*sizeof(char));
	strcpy(temp,PosString);
	length=0;
	num=split(part,temp,",");
	for(i=0;i<num;i++)
	{
		split(pos,part[i],"..");
		start=getint(pos[0]);
		end=getint(pos[1]);
		length=length+end-start+1;
	}
	return length;
}

struct GenePos
{
	char *Name;
	char *Genome;
	int rRNA_flag;
	char *Pos;
	struct GenePos *Next;
};

struct GenePos *Create_GenePos(char id[],char Genome[],int flag,char PosString[])
{
	int size,length;
	struct GenePos *node;

	size=sizeof(struct GenePos);
	node=(struct GenePos *)malloc(size);
	length=strlen(id);
	length++;
	node->Name=(char *)malloc(length*sizeof(char));
	memset(node->Name,'\0',length*sizeof(char));
	strcpy(node->Name,id);

	length=strlen(Genome);
	length++;
	node->Genome=(char *)malloc(length*sizeof(char));
	memset(node->Genome,'\0',length*sizeof(char));
	strcpy(node->Genome,Genome);

	length=strlen(PosString);
	length++;
	node->Pos=(char *)malloc(length*sizeof(char));
	memset(node->Pos,'\0',length*sizeof(char));
	strcpy(node->Pos,PosString);

	node->rRNA_flag=flag;
	node->Next=NULL;
	
	return node;
}	
	
struct GenePos *ReadGenbank(char *file,int min)
{
	FILE *fp;
	char readIn[200],Genome[100],*gene[10],id[100],PosString[2000];
	struct GenePos *head,*node,*new;
	int length,i,flag,num,yes,pos_flag,FirstFlag;

	fp=fopen(file,"r");
	if(fp==NULL)
	{
		printf("Can't open %s file!\n",file);
		exit(1);
	}
	flag=0;
	num=0;
	FirstFlag=0;
	yes=0;
	pos_flag=0;
	for(i=0;i<10;i++)
	{
		gene[i]=malloc(100*sizeof(char));
		memset(gene[i],'\0',100*sizeof(char));
	}
	memset(readIn,'\0',200*sizeof(char));
	while(fgets(readIn,200*sizeof(char),fp)!=NULL)
	{
		if(memcmp(readIn,"VERSION",7)==0)
		{
			length=strlen(readIn);
			readIn[length-1]='\0';
			i=7;
			while(readIn[i]==' ')
				i++;
			memset(Genome,'\0',100*sizeof(char));
			strcpy(Genome,&readIn[i]);
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
		if(checkFlag(readIn,"gene  ",1)==1)
		{
			flag=1;
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
		if(flag&&(checkFlag(readIn,"/pseudo",0)==1))
		{
			flag=0;
			num--;
			memset(gene[num],'\0',100*sizeof(char));
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
	//gene name
		if(flag==1&&checkFlag(readIn,"/gene=",0)==1)
		{
			memset(id,'\0',100*sizeof(char));
			length=strlen(readIn);
			readIn[length-2]='\0';
			i=0;
			while(readIn[i]!='"')
				i++;
			strcpy(id,&readIn[i+1]);

			if(num==10) //exceed
			{
				for(i=0;i<9;i++)
				{
					memset(gene[i],'\0',100*sizeof(char));
					strcpy(gene[i],gene[i+1]);
				}
				num--;
				memset(gene[num],'\0',100*sizeof(char));
			}
			strcpy(gene[num],id);
			flag++;
			num++;
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
		if(flag==1&&checkFlag(readIn,"/locus_tag=",0)==1)
		{
			memset(id,'\0',100*sizeof(char));
			length=strlen(readIn);
			readIn[length-2]='\0';
			i=0;
			while(readIn[i]!='"')
				i++;
			strcpy(id,&readIn[i+1]);

			if(num==10)
			{
				for(i=0;i<9;i++)
				{
					memset(gene[i],'\0',100*sizeof(char));
					strcpy(gene[i],gene[i+1]);
				}
				num--;
				memset(gene[num],'\0',100*sizeof(char));
			}
			strcpy(gene[num],id);
			flag++;
			num++;
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
	//gene is mRNA, CDS or ncRNA
		if(checkFlag(readIn,"mRNA  ",1)==1||checkFlag(readIn,"ncRNA  ",1)==1||checkFlag(readIn,"CDS  ",1)==1)
		{
			yes=1;
			i=0;
			while(readIn[i]==' ')
				i++;
			while(readIn[i]!=' ')
				i++;
			while(readIn[i]==' ')
				i++;
			length=strlen(readIn);
			readIn[length-1]='\0';
			memset(PosString,'\0',2000*sizeof(char));
			strcpy(PosString,&readIn[i]);
			pos_flag=1;
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
	//pos string
		if(pos_flag)
		{
			if(readIn[0]!=' ')
			{
				pos_flag=0;
				memset(readIn,'\0',200*sizeof(char));
				continue;
			}
			i=0;
			while(readIn[i]==' ')
				i++;
			if(readIn[i]=='/')
				pos_flag=0;
			else
			{
				length=strlen(readIn);
				readIn[length-1]='\0';
				strcat(PosString,&readIn[i]);
				memset(readIn,'\0',200*sizeof(char));
				continue;
			}
		}
		if(yes&&(checkFlag(readIn,"/gene=",0)==1))
		{
			if(num==0)  //wrong mRNA
			{
				memset(readIn,'\0',200*sizeof(char));
				yes=0;
				continue;
			}
			memset(id,'\0',100*sizeof(char));
			length=strlen(readIn);
			readIn[length-2]='\0';
			i=0;
			while(readIn[i]!='"')
				i++;
			strcpy(id,&readIn[i+1]);
			if(strcmp(gene[num-1],id)==0)
			{
				num--;
				memset(gene[num],'\0',100*sizeof(char));
				if(getLength(PosString)<min)
				{
					memset(readIn,'\0',200*sizeof(char));
					yes=0;
					continue;
				}
				new=Create_GenePos(id,Genome,(yes-1),PosString);
				if(FirstFlag==0)
				{
					FirstFlag++;
					head=new;
					node=new;
				}
				else
				{
					node->Next=new;
					node=new;
				}
			}
			else
			{
				for(i=0;i<num;i++)
				{
					if(strcmp(gene[i],id)==0)
					{
						if(getLength(PosString)<min)
						{
							yes=0;
							break;
						}
						new=Create_GenePos(id,Genome,(yes-1),PosString);
						if(FirstFlag==0)
						{
							FirstFlag++;
							head=new;
							node=new;
						}
						else
						{
							node->Next=new;
							node=new;
						}
						memset(readIn,'\0',200*sizeof(char));
						break;
					}
				}
			//if i=num: wrong mRNA
				if(i!=num)
				{
					for(;i<num-1;i++)
					{
						memset(gene[i],'\0',100*sizeof(char));
						strcpy(gene[i],gene[i+1]);
					}
					num--;
					memset(gene[num],'\0',100*sizeof(char));
				}
			}
			memset(readIn,'\0',200*sizeof(char));
			yes=0;
			continue;
		}
		if(yes&&checkFlag(readIn,"/locus_tag=",0)==1)
		{
			if(num==0)  //wrong mRNA
			{
				memset(readIn,'\0',200*sizeof(char));
				yes=0;
				continue;
			}
			memset(id,'\0',100*sizeof(char));
			length=strlen(readIn);
			readIn[length-2]='\0';
			i=0;
			while(readIn[i]!='"')
				i++;
			strcpy(id,&readIn[i+1]);
			if(strcmp(gene[num-1],id)==0)
			{
				num--;
				memset(gene[num],'\0',100*sizeof(char));
				if(getLength(PosString)<min)
				{
					memset(readIn,'\0',200*sizeof(char));
					yes=0;
					continue;
				}
				new=Create_GenePos(id,Genome,(yes-1),PosString);
				if(FirstFlag==0)
				{
					FirstFlag++;
					head=new;
					node=new;
				}
				else
				{
					node->Next=new;
					node=new;
				}
			}
			else
			{
				for(i=0;i<num;i++)
				{
					if(strcmp(gene[i],id)==0)
					{
						if(getLength(PosString)<min)
						{
							yes=0;
							break;
						}
						new=Create_GenePos(id,Genome,(yes-1),PosString);
						if(FirstFlag==0)
						{
							FirstFlag++;
							head=new;
							node=new;
						}
						else
						{
							node->Next=new;
							node=new;
						}
						break;
					}
				}
			//if i!=$num: wrong mRNA
				if(i!=num)
				{
					for(;i<num-1;i++)
					{
						memset(gene[i],'\0',100*sizeof(char));
						strcpy(gene[i],gene[i+1]);
					}
					num--;
					memset(gene[num],'\0',100*sizeof(char));
				}
			}
			memset(readIn,'\0',200*sizeof(char));
			yes=0;
			continue;
		}
	//rRNA
		if(checkFlag(readIn,"rRNA",1)==1)
		{
			flag=0;
			yes=2;
			i=0;
			while(readIn[i]==' ')
				i++;
			while(readIn[i]!=' ')
				i++;
			while(readIn[i]==' ')
				i++;
			length=strlen(readIn);
			readIn[length-1]='\0';
			memset(PosString,'\0',2000*sizeof(char));
			strcpy(PosString,&readIn[i]);
			if(getLength(PosString)<min)
			{
				memset(readIn,'\0',200*sizeof(char));
				yes=0;
				continue;
			}
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
	//ensure the right $id
		if(readIn[0]!=' ')
		{
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
		i=0;
		length=strlen(readIn);
		while((readIn[i]==' ')&&(i<length))
			i++;
		while((readIn[i]!=' ')&&(i<length))
			i++;
		if(i>=length)
		{
			memset(readIn,'\0',200*sizeof(char));
			continue;
		}
		while(i<length-1)
		{
			if(readIn[i]!='.')
			{
				i++;
				continue;
			}
			if(readIn[i+1]!='.')
			{
				i++;
				continue;
			}
			if((readIn[i-1]<'0')||(readIn[i-1]>'9'))
			{
				i++;
				continue;
			}
			flag=0;
			yes=0;
			break;
		}
		memset(readIn,'\0',200*sizeof(char));
	}
	fclose(fp);
	for(i=0;i<10;i++)
		free(gene[i]);
	return head;
}

struct GenePos *ReadGff(char *file, int min)
{
	FILE *fp;
	char readIn[4000],*part[10],Genome[100],id[100],PosString[2000];
	int flag,length,CDS_flag,RNA_flag,num,i,FirstFlag,j,match_flag;
	struct GenePos *head,*new,*node;

	fp=fopen(file,"r");
	if(fp==NULL)
	{
		printf("Can't open %s file!\n",file);
		exit(1);
	}
	length=0;
	flag=0;
	CDS_flag=0;
	RNA_flag=0;
	FirstFlag=0;
	match_flag=0;
	memset(readIn,'\0',4000*sizeof(char));
	while(fgets(readIn,4000*sizeof(char),fp)!=NULL)
	{
//printf("%s\n",readIn);
		if(readIn[0]=='#')
		{
			memset(readIn,'\0',4000*sizeof(char));
			continue;
		}
		num=split(part,readIn,"\t");

	//some line may be empty
		if(num<9)
		{
			memset(readIn,'\0',4000*sizeof(char));
			continue;
		}
	//printf("A0\n");
		if(strcmp(part[2],"gene")==0)
		{
			CDS_flag=0;
			RNA_flag=0;
			if(length>=min)
			{
				i=strlen(PosString);
				if(PosString[i-1]==',')
					PosString[i-1]='\0';
				new=Create_GenePos(id,Genome,0,PosString);
				if(FirstFlag==0)
				{
					FirstFlag++;
					head=new;
					node=new;
				}
				else
				{
					node->Next=new;
					node=new;
				}
			}
		//gene name
			i=0;
			while(memcmp(&part[8][i],"Name=",5)!=0)
				i++;
			j=i+5;
			while(part[8][j]!=';'&&part[8][j]!='\n')
				j++;
			part[8][j]='\0';
			memset(id,'\0',100*sizeof(char));
			strcpy(id,&part[8][i+5]);
			memset(Genome,'\0',100*sizeof(char));
			strcpy(Genome,part[0]);
			memset(PosString,'\0',2000*sizeof(char));
			length=0;
			flag=1;
			memset(readIn,'\0',4000*sizeof(char));
			continue;
		}
	//printf("A1\n");
		if(strcmp(part[2],"pseudogene")==0)
		{
			flag=2;
			CDS_flag=0;
			RNA_flag=0;
			memset(readIn,'\0',4000*sizeof(char));
			continue;
		}
	//printf("A2\n");
		if(flag==1)
		{
			match_flag=1;
			j=strlen(part[8]);
			i=0;
			while(i<j&&memcmp(&part[8][i],"gbkey=",6)!=0)
				i++;
			if(i==j)
			{
				CDS_flag=0;
				RNA_flag=0;
				match_flag=0;
			}
			i=i+6;
			if(match_flag&&memcmp(&part[8][i],"mRNA",4)!=0&&memcmp(&part[8][i],"ncRNA",5)!=0&&memcmp(&part[8][i],"CDS",3)!=0)
			{
				CDS_flag=0;
				RNA_flag=0;
				match_flag=0;
			}
			if(match_flag)
			{
				flag++;
				length=0;
				memset(PosString,'\0',2000*sizeof(char));
				if(part[6][0]=='-')
					strcat(PosString,"complement:");
				if(memcmp(&part[8][i],"CDS",3)==0)
				{
					CDS_flag=1;
					length=atoi(part[4])-atoi(part[3])+1;
					strcat(PosString,part[3]);
					strcat(PosString,"..");
					strcat(PosString,part[4]);
					strcat(PosString,",");
				}
				else
					RNA_flag=1;
				memset(readIn,'\0',4000*sizeof(char));
				continue;
			}//else: go on
		}
	//printf("A3\n");
		if(strcmp(part[2],"exon")==0)
		{
			if(RNA_flag)
			{
				length=length+atoi(part[4])-atoi(part[3])+1;
				strcat(PosString,part[3]);
				strcat(PosString,"..");
				strcat(PosString,part[4]);
				strcat(PosString,",");
				memset(readIn,'\0',4000*sizeof(char));
				continue;
			}
			memset(readIn,'\0',4000*sizeof(char));
			continue;
		}
	//printf("A4\n");
		if(strcmp(part[2],"CDS")==0)
		{
			if(CDS_flag)
			{
				length=length+atoi(part[4])-atoi(part[3])+1;
				strcat(PosString,part[3]);
				strcat(PosString,"..");
				strcat(PosString,part[4]);
				strcat(PosString,",");
				memset(readIn,'\0',4000*sizeof(char));
				continue;
			}
			memset(readIn,'\0',4000*sizeof(char));
			continue;
		}
	//printf("A5\n");
	//rRNA
		if(strcmp(part[2],"rRNA")==0)
		{
			if(length>=min)
			{
				i=strlen(PosString);
				if(PosString[i-1]==',')
					PosString[i-1]='\0';
				new=Create_GenePos(id,Genome,0,PosString);
				if(FirstFlag==0)
				{
					FirstFlag++;
					head=new;
					node=new;
				}
				else
				{
					node->Next=new;
					node=new;
				}
			}
			flag=2;
			length=atoi(part[4])-atoi(part[3])+1;
			if(length>=min)
			{
				memset(Genome,'\0',100*sizeof(char));
				strcpy(Genome,part[0]);
				memset(PosString,'\0',2000*sizeof(char));
				if(part[6][0]=='-')
					strcat(PosString,"complement:");
				strcat(PosString,part[3]);
				strcat(PosString,"..");
				strcat(PosString,part[4]);
				new=Create_GenePos(id,Genome,1,PosString);
				if(FirstFlag==0)
				{
					FirstFlag++;
					head=new;
					node=new;
				}
				else
				{
					node->Next=new;
					node=new;
				}
			}
			length=0;
			RNA_flag=0;
			CDS_flag=0;
			memset(readIn,'\0',4000*sizeof(char));
			continue;
		}
	//printf("A6\n");
	//reset the flag
		RNA_flag=0;
		CDS_flag=0;
		memset(readIn,'\0',4000*sizeof(char));
	}
	fclose(fp);
	if(length>=min)
	{
		i=strlen(PosString);
		if(PosString[i-1]==',')
			PosString[i-1]='\0';
		new=Create_GenePos(id,Genome,0,PosString);
		if(FirstFlag==0)
		{
			FirstFlag++;
			head=new;
			node=new;
		}
		else
		{
			node->Next=new;
			node=new;
		}
	}
	return head;
}

int GetGenomeSequence(FILE *fp,char *Genome,char *Genome_name)
{
	char read[300];
	int flag,success,length;
	
	flag=0;
	success=0;
	memset(read,'\0',300*sizeof(char));
	length=strlen(Genome_name);
	while(fgets(read,300*sizeof(char),fp)!=NULL)
	{
		if(read[0]=='>')
		{
			if(memcmp(&read[1],Genome_name,length)==0)
			{
				flag=1;
				break;
			}
		}
		memset(read,'\0',300*sizeof(char));
	}
	if(flag==0)
		return 0; //fail

	while(fgets(read,300*sizeof(char),fp)!=NULL)
	{
		if(read[0]=='>')
		{
			fseek(fp,-10,1);
			break;
		}
		length=strlen(read);
		read[length-1]='\0';
		strcat(Genome,read);
	}
	return 1;
}

int GeneratecDNA(char *Genome,char *Pos,char *cDNA,int cDNA_length)
{
	int length,start,end,i,j,num;
	char temp[2000],*cDNA_temp; //maybe don't need, it's the last time to use Pos
	char *part[200],*pos[2];

	memset(temp,'\0',2000*sizeof(char));
	strcpy(temp,Pos);
	length=strlen(Genome);

	cDNA_temp=(char *)malloc((cDNA_length+1)*sizeof(char));
	memset(cDNA_temp,'\0',(cDNA_length+1)*sizeof(char));
	j=0;
	num=split(part,temp,",");
	for(i=0;i<num;i++)
	{
		split(pos,part[i],"..");
		start=getint(pos[0]);
		end=getint(pos[1]);
		if(end>length)
			return 0;

		while(start<=end)
		{
			cDNA_temp[j]=Genome[start-1];
			j++;
			start++;
		}		
	}
	if(strstr(Pos,"complement")!=NULL)
		minus_strand(cDNA_temp,cDNA,cDNA_length);
	else
	{
		for(i=0;i<cDNA_length;i++)
			cDNA[i]=cDNA_temp[i];
	}
	free(cDNA_temp);
	return 1;
}

void ReadGeneList(char *file,struct Genome *headG)
{
	FILE *fp;
	char readIn[400];
	int FirstFlag,num,length,size;
	char *part[10];
	struct Genome *nodeG;
	struct InputGene *head,*new,*node;

	fp=fopen(file,"r");
	if(fp==NULL)
	{
		printf("Error! Can't open %s file!\n",file);
		exit(1);
	}
	FirstFlag=0;
	size=sizeof(struct InputGene);
	memset(readIn,'\0',400*sizeof(char));
	while(fgets(readIn,400*sizeof(char),fp)!=NULL)
	{
		if(readIn[0]=='[')
		{
			if(FirstFlag==2)
				nodeG->gene_list=head;
			length=strlen(readIn);
			readIn[length-2]='\0';
			nodeG=headG;
			FirstFlag=0;
			while(nodeG)
			{
				if(strcmp(nodeG->Name,&readIn[1])==0)
				{
					FirstFlag=1;
					break;
				}
				nodeG=nodeG->Next;
			}
			if(FirstFlag==0)
				printf("Warning: the %s genome in %s file is not in the abundance table!\n",&readIn[1],file);
			memset(readIn,'\0',400*sizeof(char));
			continue;
		}
		if(FirstFlag==0)
		{
			memset(readIn,'\0',400*sizeof(char));
			continue;
		}
		length=strlen(readIn);
		readIn[length-1]='\0';
		num=split(part,readIn,"\t");

		if(num!=2&&num!=3)
		{
			memset(readIn,'\0',400*sizeof(char));
			continue;
		}
		
		new=(struct InputGene *)malloc(size);
		memset(new->Name,'\0',100*sizeof(char));
		strcpy(new->Name,part[0]);
		new->count=atoi(part[1]);
		if(num==3)
			new->FC=atof(part[2]);
		else
			new->FC=1;
		new->Next=NULL;
		if(FirstFlag==1)
		{
			head=new;
			node=new;
			FirstFlag++;
		}
		else
		{
			node->Next=new;
			node=new;
		}
	}
	fclose(fp);
	if(FirstFlag==2)
		nodeG->gene_list=head;
}

int GenerateRandomCount(float probability[])
{
	float random;
	int i;

	random=1.0*rand()/RAND_MAX;
	for(i=1;i<=1000;i++)
		if(random<=probability[i])
			return i;
	return 1;
}

struct DNASequences
{
	char *Name;
	char *Fasta;
	struct DNASequences *Next;
};

struct DNASequences *ReadFastaFile(char *file,int min)
{
	int size,i,j,k,FirstFlag,length;
	char *readIn;
	FILE *fp;
	struct DNASequences *head,*new,*node;

	fp=fopen(file,"r");
	if(fp==NULL)
	{
		printf("Error! Can't open %s file!\n",file);
		return NULL;
	}

	size=FileSize(file);
	size++;
	readIn=(char *)malloc(size*sizeof(char));
	memset(readIn,'\0',size*sizeof(char));
	fread(readIn,size*sizeof(char),1,fp);
	fclose(fp);

	i=0;
	FirstFlag=0;
	while(i<size)
	{
		if(readIn[i]=='>')
		{
		//name
			j=i+1;
			while(j<size&&readIn[j]!=' '&&readIn[j]!='\n')
				j++;
			if(j==size)
				return head;
			if(readIn[j]=='\n')
				k=j;  //start sequence
			else
			{
				k=j+1;
				while(k<size&&readIn[k]!='\n')
					k++;
			}
			readIn[j]='\0';
			new=(struct DNASequences *)malloc(sizeof(struct DNASequences));
			new->Name=(char *)malloc((j-i)*sizeof(char));
			memset(new->Name,'\0',(j-i)*sizeof(char));
			strcpy(new->Name,&readIn[i+1]);

		//sequence
			i=k+1;
			j=i;
			length=0;
			while(j<size&&readIn[j]!='>')
			{
				if((readIn[j]!=' ')&&(readIn[j]!='\n'))
					length++;
				j++;
			}
			if(length<min)
			{
				free(new->Name);
				free(new);
				i=j;
				continue;
			}
			new->Fasta=(char *)malloc((j-i+1)*sizeof(char));
			memset(new->Fasta,'\0',(j-i+1)*sizeof(char));
			k=i;
			while(k<j)
			{
				if(readIn[k]=='\n')
				{
					readIn[k]='\0';
					strcat(new->Fasta,&readIn[i]);
					i=k+1;
				}
				k++;
			}
			i=j;
			new->Next=NULL;
			if(FirstFlag==0)
			{
				head=new;
				node=new;
				FirstFlag++;
			}
			else
			{
				node->Next=new;
				node=new;
			}
		}
		else
			i++;
	}
	free(readIn);
	return head;
}

void usage()
{
	printf("This is a customizable metatranscriptome simulation system: NeSSMt (Next-generation Sequencing Simulator for Metatranscriptomics).\n\n");
	printf("Usage: ./NeSSMt -abundance abundance_file -index index_file -o prefix_output [options]\n");
	printf("[Necessary parameter]:\n");
	printf("-abundance <abundance_file>     [char]this file contains the organisms' name and their abundance.\n");
	printf("-index <index_file>             [char]the index file contains the information of reference genomes, transcripts and so on, it is generated by mk_index.pl script.\n");
	printf("-o <output_file>                [char]name of output file.\n\n");

	printf("[Optional parameter]:\n");
	printf("-r <reads_number>               [int]total number of simulated reads, default 10000.\n");
	printf("-single                         simulate single reads. By default NeSSMt simulates paired reads.\n");
	printf("-l <read_length>                [int]length of the reads, default 150bp.\n");
	printf("-fragment <fragment_length>     [int]length of fragment in paired-end sequencing, default 350bp.\n");
	printf("-sd <standard deviation>        [int]standard deviation of fragment length, default 20.\n");
	printf("-config <config_file>           [char]the configure file for sequencing errors, default is the \"simulation.config\" file.\n");
	printf("-seed <seed>                    [int]the seed used in random function, default 0.\n");
	printf("-strandspecific                 simulate strand-specific RNA-seq sequencing, default is unstrand-specific.\n");
	printf("-express <express_per>          [float]the percentage of expressed gene in all genes for each organism, default 0.1.\n");
	printf("-gene <gene_file>               [char]this file contains the list of expressed genes and their relative abundance.\n");
	printf("-rRNA <percentage>		[float]the percentage of reads from rRNA genes in all simulated reads, default 0.1.\n\n");

	printf("[Simulate differential expression]:\n");
	printf("-diff				simulate differential expression.\n");
	printf("-up <up_expression_per>         [float]the percentage of up regulated genes in all expressed genes for each organism, default 0.1.\n");
	printf("-down <down_expression_per>     [float]the percentage of down regulated genes in all expressed genes for each organism, default 0.1.\n");
	printf("-maxFC <max_fold change>        [float]the max fold change in differential expression(>1), default 20.\n");
	printf("-minFC <min_fold change>        [float]the min fold change in differential expression(>1), default 2.\n");
	printf("-more <more_abundance_per>	[float]the percentage of organisms whose relative abundance are more in control group than in treat group, default 0.2.\n");
	printf("-less <down_abundance_per>	[float]the percentage of organisms whose relative abundance are less in control group than in treat group, default 0.2.\n");
	printf("-maxAd <max_abundance_change>	[float]the max fold change in relative abundance(>1), default 10.\n");
	printf("-minAd <min_abundance_change>	[float]the min fold change in relative abundance(>1), default 2.\n\n");

	printf("[Simulate biological replicates]:\n");
	printf("-replicate <replicates>         [int]the number of biological replicates. By default, replicates are 3 in differential expression simulation and replicates is 1 for other simulation.\n");
	printf("-prob <probability>             [float]the probability is used to simulate expression counts in biological replicates, the bigger of probability, the more dispersion of counts, default 0.33333.\n");
	return;
}

int check_input(char abundance_file[],char index_file[],char output_prefix[],int readNum,int Len_fragment,int Len_read,char config_file[],float express_per,char gene_file[],float up_expression_per,float down_expression_per,float minFC,float maxFC,float pnbinom,float rRNA,float more_abundance_per,float less_abundance_per,float maxAbundance,float minAbundance)
{
	if(strlen(abundance_file)==0)
	{
		printf("Error! Please input the abundance file by \"-abundance\" parameter.\n");
		usage();
		return 0;
	}
	if(access(abundance_file,F_OK)!=0)
	{
		printf("Error! Can't find the abundance file: %s.\n",abundance_file);
		return 0;
	}
	if(strlen(index_file)==0)
	{
		printf("Error! Please input the index file by \"-index\" parameter.\n");
		usage();
		return 0;
	}
	if(access(index_file,F_OK)!=0)
	{
		printf("Error! Can't find the index file: %s.\n",index_file);
		return 0;
	}
	if(strlen(output_prefix)==0)
	{
		printf("Error! the prefix of output file is needed, you can input this prefix by \"-o\" parameter.\n");
		usage();
		return 0;
	}
	if(readNum<=0)
	{
		printf("Error! the total number of simulated reads must be larger than 0, now is %d.\n",readNum);
		printf("  You can fix it by \"-r\" parameter.\n");
		return 0;
	}
	if(Len_fragment<=Len_read)
	{
		printf("Error! the length of fragment must be bigger than the length of sequencing read.\n");
		printf("  Now the length of fragment is %d bp, read length is %d bp.\n",Len_fragment,Len_read);
		usage();
		return 0;
	}
	if(access(config_file,F_OK)!=0)
	{
		printf("Error! Can't find the %s file!\n",config_file);
		printf("  You can fix it by \"-config\" parameter.\n");
		usage();
		return 0;
	}
	if(express_per>=1||express_per<=0)
	{
		printf("Error! The percentage of expression gene in one organism should be between 0 and 1, now is %f.\n",express_per);
		printf("  You can fix it by \"-express\" parameter.\n");
		return 0;
	}
	if(strlen(gene_file)>0&&(access(gene_file,F_OK)!=0))
	{
		printf("Error! Can't find the %s file.\n",gene_file);
		printf("  You can fix it by \"-gene\" parameter.\n");
		usage();
		return 0;
	}
	if(rRNA>1||rRNA<0)
	{
		printf("Error! The percentage of rRNA reads in all simulated reads should be between 0 and 1, now is %f.\n",rRNA);
		printf("  You can fix it by \"-rRNA\" parameter.\n");
		return 0;
	}
	if(up_expression_per>=1||up_expression_per<=0)
	{
		printf("Error! The percentage of up expression genes in all expressed genes should be between 0 and 1, now is %f.\n",up_expression_per);
		printf("  You can fix it by \"-up\" parameter.\n");
		return 0;
	}
	if(down_expression_per>=1||down_expression_per<=0)
	{
		printf("Error! The percentage of down expression genes in all expressed genes should be between 0 and 1, now is %f.\n",down_expression_per);
		printf("  You can fix it by \"-down\" parameter.\n");
		return 0;
	}
	if(up_expression_per+down_expression_per>=1)
	{
		printf("Error! The total percentage of up expression genes and down expression genes must be less than 1.\n");
		printf("  Now the total percentage is %f.\n",(up_expression_per+down_expression_per));
		return 0;
	}
	if(minFC<=1)
	{
		printf("Error! The min fold change should be larger than 1, now is %f.\n",minFC);
		printf("  You can fix it by \"-minFC\" parameter.\n");
		return 0;
	}
	if(maxFC<minFC)
	{
		printf("Error! The max fold change should be larger than the min fold change.\n");
		printf("  Now the max fold change is %f, min fold change is %f.\n",maxFC,minFC);
		return 0;
	}
	if(more_abundance_per>=1||more_abundance_per<=0)
	{
		printf("Error! The percentage of genoms whose relative abundance in control are more than in treat group should be between 0 and 1, now is %f.\n",more_abundance_per);
		printf("  You can fix it by \"-more\" parameter.\n");
		return 0;
	}
	if(less_abundance_per>=1||less_abundance_per<=0)
	{
		printf("Error! The percentage of genoms whose relative abundance in control are less than in treat group should be between 0 and 1, now is %f.\n",less_abundance_per);
		printf("  You can fix it by \"-less\" parameter.\n");
		return 0;
	}
	if(down_expression_per>=1||down_expression_per<=0)
	{
		printf("Error! The percentage of down expression genes in all expressed genes should be between 0 and 1, now is %f.\n",down_expression_per);
		printf("  You can fix it by \"-down\" parameter.\n");
		return 0;
	}
        if(more_abundance_per+less_abundance_per>=1)
	{
		printf("Error! The total percentage of genomes whose relative abundance changed must be less than 1.\n");
		printf("  Now the total percentage is %f.\n",(more_abundance_per+less_abundance_per));
		return 0;
	}
	if(minAbundance<=1)
	{
		printf("Error! The min fold change in relative abundance should be larger than 1, now is %f.\n",minAbundance);
		printf("  You can fix it by \"-minAd\" parameter.\n");
		return 0;
	}
	if(maxAbundance<minAbundance)
	{
		printf("Error! The max fold change in relative abundance should be larger than the min fold change.\n");
		printf("  Now the max fold change in relative abundance is %f, min fold change is %f.\n",maxAbundance,minAbundance);
		return 0;
	}
	if(pnbinom>1||pnbinom<=0)
	{
		printf("Error! The probability used in the negative binomial should be between 0 and 1, now is %f.\n",pnbinom);
		printf("  You can fix it by \"-prob\" parameter.\n");
		return 0;
	}
	return 1;
}

int Create_output(FILE *fp[],int turn,char *output_prefix,char *Control_treat,int Pair_turn)
{
	int length;
	char *output,string[5];

	length=strlen(output_prefix);
	length=length+20;
	output=(char *)malloc(length*sizeof(char));
	strcpy(output,output_prefix);
	if(Control_treat!=NULL)
		strcat(output,Control_treat);
	if(turn!=0)
	{
		sprintf(string,"%d",turn);
		strcat(output,string);
	}
	if(Pair_turn==0)
		strcat(output,".fq");
	else if(Pair_turn==1)
		strcat(output,"_1.fq");
	else
		strcat(output,"_2.fq");

	if(turn==0)
	{
		fp[turn]=fopen(output,"w");
		if(fp[turn]==NULL)
		{
			printf("Error! Can't create %s file.\n",output);
			free(output);
			return 0;
		}
	}
	else
	{
		fp[turn-1]=fopen(output,"w");
		if(fp[turn-1]==NULL)
		{
			printf("Error! Can't create %s file.\n",output);
			free(output);
			return 0;
		}
	}
	free(output);
	return 1;
}

void CleanMemory(struct Genome *headGenome,float *err_rand,int replicates,gsl_rng *r,FILE *T1[],FILE *T2[],FILE *C1[],FILE *C2[],int Flag_control_treat,int Flag_pair)
{
	struct Genome *nodeGenome;
	struct InputGene *nodeGeneList,*tempGeneList;
	int i;

	free(err_rand);
	while(headGenome)
	{
		nodeGenome=headGenome;
		headGenome=headGenome->Next;
		if(nodeGenome->transcript_path!=NULL)
			free(nodeGenome->transcript_path);
		if(nodeGenome->genome_path!=NULL)
			free(nodeGenome->genome_path);
		if(nodeGenome->annotation_path!=NULL)
			free(nodeGenome->annotation_path);
		if(nodeGenome->rRNA_path!=NULL)
			free(nodeGenome->rRNA_path);
		if(nodeGenome->gene_list!=NULL)
		{
			nodeGeneList=nodeGenome->gene_list;
			while(nodeGeneList)
			{
				tempGeneList=nodeGeneList->Next;
				free(nodeGeneList);
				nodeGeneList=tempGeneList;
			}
		}
		free(nodeGenome);
	}

//if repeat, free
	if(replicates>1)
	{
		gsl_rng_free (r);
		for(i=1;i<=replicates;i++)
		{
			if(Flag_control_treat==1)
			{
				if(T1[i-1]!=NULL)
					fclose(T1[i-1]);
				if(C1[i-1]!=NULL)
					fclose(C1[i-1]);
				if(Flag_pair==1)
				{
					if(T2[i-1]!=NULL)
						fclose(T2[i-1]);
					if(C2[i-1]!=NULL)
						fclose(C2[i-1]);
				}
			}
			else
			{
				if(T1[i-1]!=NULL)
					fclose(T1[i-1]);
				if(Flag_pair==1&&(T2[i-1]!=NULL))
					fclose(T2[i-1]);
			}
		}
	}
	else
	{
		if(Flag_control_treat==1)
		{
			if(T1[0]!=NULL)
				fclose(T1[0]);
			if(C1[0]!=NULL)
				fclose(C1[0]);
			if(Flag_pair==1)
			{
				if(T2[0]!=NULL)
					fclose(T2[0]);
				if(C2[0]!=NULL)
					fclose(C2[0]);
			}
		}
		else
		{
			if(T1[0]!=NULL)
				fclose(T1[0]);
			if(Flag_pair==1&&(T2[0]!=NULL))
				fclose(T2[0]);
		}
	}
}

int rand_insert(int average,int sd,int max,int min)
{
	float out,random,test;
	int flag=0;

	if(sd==0)
		return average;

	while(flag==0)
	{
		random=(float)((float)rand()/RAND_MAX-0.5)*50*sd+average;
		test=(float)rand()/RAND_MAX*2.0/sqrt(2.0*3.1415926)/sd;
		out=random;
		random=1.0/sqrt(2.0*3.1415926)/sd*(float)exp((-1)*(random-average)*(random-average)/2.0/sd/sd);

		if((test<random)&&(out>min)&&(out<=max))    //round-off
			flag=1;
	}
	return (int)(out+0.5);
}

char sub_single(char b1,char b2,char b3,float p1,float p2,float p3)
{
	float per1,per2,per;
	per=(float)(rand())/(float)(RAND_MAX);
	per1=p1/(p1+p2+p3);
	per2=p2/(p1+p2+p3);
	if(per<=per1)
		return b1;
	else if(per>per1 && per<=(per1+per2))
		return b2;
	else
		return b3;
}
//======================================================================================
//function:sub_all--substitute the inputed base
//======================================================================================
char sub_all(char seq,float ac,float ag,float at,float cg,float ct,float ca,float gt,float ga,float gc,float ta,float tc,float tg)
{
	if(seq=='A')
		return sub_single('C','G','T',ac,ag,at);
	else if(seq=='C')
		return sub_single('G','T','A',cg,ct,ca);
	else if(seq=='G')
		return sub_single('T','A','C',gt,ga,gc);
	else if(seq=='T')
		return sub_single('A','C','G',ta,tc,tg);
	return 'N'; //in order to avoid warnings in gcc
}

//======================================================================================
//function:ins_all--generate a base for insertion
//======================================================================================
char ins_all()
{
	float random=(float)(rand())/(float)(RAND_MAX);
	if(random<=0.25)
		return 'A';
	else if(random<=0.5)
		return 'C';
	else if(random<=0.75)
		return 'G';
	else
		return 'T';
}

int Error_record(int position_err,int Error_point,char a1,char a2,char *Error)
{
	char s[6];
	sprintf(s,"%d",position_err);
	if(Error_point!=0)
	{
		Error[Error_point]='|';
		Error_point++;
	}
	strcat(&Error[Error_point],s);
	Error_point+=strlen(s);
	Error[Error_point]=':';
	Error_point++;
	Error[Error_point]=a1;
	Error_point++;
	Error[Error_point]=':';
	Error_point++;
	Error[Error_point]=a2;
	Error_point++;
	return Error_point;
}

void bp(float type[],char *Read,char *Phred,char *Ref,float subratio[],char *Error,int Read_length,float *err_rand,int position,float qual[],int Ref_length)
{
	int a,turn,Error_point,add_mark=0;
	float random,err_random,p;
	int quality();

	Error_point=0;
	for(a=0;(a+add_mark<Read_length)&&((a+add_mark+position)<Ref_length);a++)
	{
		random=(float)(rand())/(float)(RAND_MAX);
		turn=quality((a+add_mark),err_rand);
		p=qual[turn];
		if(random<=p)
		{
			err_random=(float)(rand())/(float)(RAND_MAX);
			if(err_random<=type[0])
			{
				Read[a+add_mark]=sub_all(Ref[a+position],subratio[0],subratio[1],subratio[2],subratio[3],subratio[4],subratio[5],subratio[6],subratio[7],subratio[8],subratio[9],subratio[10],subratio[11]);
				Error_point=Error_record(a+add_mark,Error_point,Ref[a+position],Read[a+add_mark],Error);
				Phred[a+add_mark]=(char)(turn+33);
			}
			else if(err_random<=type[1])
			{
				Read[a+add_mark]=ins_all();
				Error_point=Error_record(a+add_mark,Error_point,'-',Read[a+add_mark],Error);
				Phred[a+add_mark]=(char)(turn+33);
				add_mark++;
				if((a+add_mark)<Read_length)
				{
					Read[a+add_mark]=Ref[a+position];
					turn=quality((a+add_mark),err_rand);
					Phred[a+add_mark]=(char)(turn+33);
				}
			}
			else
			{
				Error_point=Error_record(a+add_mark,Error_point,Ref[a+position],'-',Error);
				add_mark--;
			}
		}
		else
		{
			Read[a+add_mark]=Ref[a+position];
			Phred[a+add_mark]=(char)(turn+33);
		}
	}
}

int NeSSM(int Flag_pair,int Flag_strand_specific,int Num_reads,int Len_read,int Len_fragment,float sd,char *cDNA,FILE *One,FILE *Two,char *Name_genome,char *Name_gene,int Turn_read,float type[],float subratio[],float *err_rand,float qual[])
{
	int Length_cDNA,i,point1,point2,Fragment_here;
	char *Read1,*Read2,*Phred1,*Phred2,Error1[500],Error2[500],*Ref_minus;

	Length_cDNA=strlen(cDNA);
	Ref_minus=(char *)malloc((Length_cDNA+1)*sizeof(char));
	memset(Ref_minus,'\0',(Length_cDNA+1)*sizeof(char));
	minus_strand(cDNA,Ref_minus,Length_cDNA);

	Read1=(char *)malloc((Len_read+1)*sizeof(char));
	Phred1=(char *)malloc((Len_read+1)*sizeof(char));
	if(Flag_pair)
	{
		Read2=(char *)malloc((Len_read+1)*sizeof(char));
		Phred2=(char *)malloc((Len_read+1)*sizeof(char));
	}
	if(Flag_pair)   //simulation pair-end
	{
		for(i=0;i<Num_reads;i++)
		{
			memset(Read1,'\0',(Len_read+1)*sizeof(char));
			memset(Read2,'\0',(Len_read+1)*sizeof(char));
			memset(Phred1,'\0',(Len_read+1)*sizeof(char));
			memset(Phred2,'\0',(Len_read+1)*sizeof(char));
			memset(Error1,'\0',500*sizeof(char));
			memset(Error2,'\0',500*sizeof(char));

			if(Flag_strand_specific||(float)(rand())/(float)(RAND_MAX)<=0.5)//minus
			{
				if(Length_cDNA<=2*Len_read)
				{
					point1=Length_cDNA-1;
					point2=0;
				}
				else
				{
					Fragment_here=rand_insert(Len_fragment,sd,Length_cDNA,Len_read);
					point1=(int)(rand()/(float)(RAND_MAX)*(Length_cDNA-Fragment_here)+Fragment_here);
					if(point1>=Length_cDNA)
						point1=Length_cDNA-1;
					point2=point1-Fragment_here+1;
					if(point2<0)
						point2=0;
				}
				bp(type,Read1,Phred1,Ref_minus,subratio,Error1,Len_read,err_rand,(Length_cDNA-point1-1),qual,Length_cDNA);
				bp(type,Read2,Phred2,cDNA,subratio,Error2,Len_read,err_rand,point2,qual,Length_cDNA);
				if(Error1[0]=='\0')
					fprintf(One,"@%d:Genome: %s|Gene: %s/1\n%s\n+pos=%d, strand=minus\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Phred1);
				else
					fprintf(One,"@%d:Genome: %s|Gene: %s/1\n%s\n+pos=%d, strand=minus, error=%s\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Error1,Phred1);
				if(Error2[0]=='\0')
					fprintf(Two,"@%d:Genome: %s|Gene: %s/2\n%s\n+pos=%d, strand=plus\n%s\n",Turn_read,Name_genome,Name_gene,Read2,point2,Phred2);
				else
					fprintf(Two,"@%d:Genome: %s|Gene: %s/2\n%s\n+pos=%d, strand=plus, error=%s\n%s\n",Turn_read,Name_genome,Name_gene,Read2,point2,Error2,Phred2);
				Turn_read++;
			}
			else //plus
			{
				if(Length_cDNA<=2*Len_read)
				{
					point1=0;
					point2=Length_cDNA-1;
				}
				else
				{
					Fragment_here=rand_insert(Len_fragment,sd,Length_cDNA,Len_read);
					point1=(int)(rand()/(float)(RAND_MAX)*(Length_cDNA-Fragment_here));
					point2=point1+Fragment_here-1;
					if(point2>=Length_cDNA)
						point2=Length_cDNA-1;
				}
				bp(type,Read1,Phred1,cDNA,subratio,Error1,Len_read,err_rand,point1,qual,Length_cDNA);
				bp(type,Read2,Phred2,Ref_minus,subratio,Error2,Len_read,err_rand,(Length_cDNA-point2-1),qual,Length_cDNA);
				if(Error1[0]=='\0')
					fprintf(One,"@%d:Genome: %s|Gene: %s/1\n%s\n+pos=%d, strand=plus\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Phred1);
				else
					fprintf(One,"@%d:Genome: %s|Gene: %s/1\n%s\n+pos=%d, strand=plus, error=%s\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Error1,Phred1);
				if(Error2[0]=='\0')
					fprintf(Two,"@%d:Genome: %s|Gene: %s/2\n%s\n+pos=%d, strand=minus\n%s\n",Turn_read,Name_genome,Name_gene,Read2,point2,Phred2);
				else
					fprintf(Two,"@%d:Genome: %s|Gene: %s/2\n%s\n+pos=%d, strand=minus, error=%s\n%s\n",Turn_read,Name_genome,Name_gene,Read2,point2,Error2,Phred2);
				Turn_read++;
			}
		}
	}
	else //simulation   single
	{
		for(i=0;i<Num_reads;i++)
		{
			memset(Read1,'\0',(Len_read+1)*sizeof(char));
			memset(Phred1,'\0',(Len_read+1)*sizeof(char));
			memset(Error1,'\0',500*sizeof(char));

			if(Flag_strand_specific||(float)(rand())/(float)(RAND_MAX)<=0.5)//minus
			{
				point1=(int)(rand()/(float)(RAND_MAX)*(Length_cDNA-Len_read)+Len_read);
				if(point1>=Length_cDNA)
					point1=Length_cDNA-1;
				bp(type,Read1,Phred1,Ref_minus,subratio,Error1,Len_read,err_rand,(Length_cDNA-point1-1),qual,Length_cDNA);
				if(Error1[0]=='\0')
					fprintf(One,"@%d:Genome: %s|Gene: %s\n%s\n+pos=%d, strand=minus\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Phred1);
				else
					fprintf(One,"@%d:Genome: %s|Gene: %s\n%s\n+pos=%d, strand=minus, error=%s\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Error1,Phred1);
				Turn_read++;
			}
			else //plus
			{
				point1=(int)(rand()/(float)(RAND_MAX)*(Length_cDNA-Len_read));
				bp(type,Read1,Phred1,cDNA,subratio,Error1,Len_read,err_rand,point1,qual,Length_cDNA);
				if(Error1[0]=='\0')
					fprintf(One,"@%d:Genome: %s|Gene: %s\n%s\n+pos=%d, strand=plus\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Phred1);
				else
					fprintf(One,"@%d:Genome: %s|Gene: %s\n%s\n+pos=%d, strand=plus, error=%s\n%s\n",Turn_read,Name_genome,Name_gene,Read1,point1,Error1,Phred1);
				Turn_read++;
			}
		}
	}
	free(Read1);
	free(Phred1);
	if(Flag_pair)
	{
		free(Read2);
		free(Phred2);
	}
	return Turn_read;
}

struct StoreGene *StoreGeneList(int Flag_control_treat,int replicates,int Expression,float FC,char *Fasta,char *Name_genome,char *Name_gene)
{
	int size;
	struct StoreGene *node;

	size=sizeof(struct StoreGene);
	node=(struct StoreGene *)malloc(size);

	node->Fasta=Fasta;
	node->Genome=Name_genome;
	node->Gene=Name_gene;
	node->Next=NULL;

	node->countTreat=(int *)malloc(replicates*sizeof(int));
	node->countTreat[0]=Expression;

	if(Flag_control_treat==1)
	{
		node->FC=FC;
		node->countControl=(int *)malloc(replicates*sizeof(int));
		if(FC==1)
			node->countControl[0]=Expression;
		else if(FC>1)
			node->countControl[0]=1.0*Expression/FC+0.5;
		else
			node->countControl[0]=-1.0*Expression*FC+0.5;
	}
	return node;
}

void Cal_replicates(int Flag_control_treat,int replicates,float probability,struct StoreGene *head,gsl_rng *r)
{
	struct StoreGene *node;
	int i,success;
	
	node=head;
	while(node)
	{
		success=1.0*node->countTreat[0]*probability+0.5;
		for(i=0;i<replicates;i++)
		{
			if(success==0)
				node->countTreat[i]=node->countTreat[0];
			else
				node->countTreat[i]=success+gsl_ran_negative_binomial(r,probability,success);
		}
		if(Flag_control_treat==1)
		{
			success=1.0*node->countControl[0]*probability+0.5;
			for(i=0;i<replicates;i++)
			{
				if(success==0)
					node->countControl[i]=node->countControl[0];
				else
					node->countControl[i]=success+gsl_ran_negative_binomial(r,probability,success);
			}
		}
		node=node->Next;
	}
	return;
}

void Simulate(int Flag_control_treat,int replicates,int Turn_Treat[],int Turn_Control[],struct StoreGene *node,int Flag_pair,int Flag_strand_specific,int Len_read,int Len_fragment,int sd,FILE *T1[],FILE *T2[],FILE *C1[],FILE *C2[],float type[],float subratio[],float *err_rand,float qual[])
{
	int i;

	for(i=0;i<replicates;i++)
		Turn_Treat[i]=NeSSM(Flag_pair,Flag_strand_specific,node->countTreat[i],Len_read,Len_fragment,sd,node->Fasta,T1[i],T2[i],node->Genome,node->Gene,Turn_Treat[i],type,subratio,err_rand,qual);
	if(Flag_control_treat==1)
	{
		for(i=0;i<replicates;i++)
			Turn_Control[i]=NeSSM(Flag_pair,Flag_strand_specific,node->countControl[i],Len_read,Len_fragment,sd,node->Fasta,C1[i],C2[i],node->Genome,node->Gene,Turn_Control[i],type,subratio,err_rand,qual);
	}
}

void Normalize(struct StoreGene *head,int replicates,int Flag_control_treat,int ReadNum)
{
	struct StoreGene *node;
	int *total,*left,i;

	total=(int *)malloc(replicates*sizeof(int));
	left=(int *)malloc(replicates*sizeof(int));

	for(i=0;i<replicates;i++)
		total[i]=0;
	node=head;
	while(node)
	{
		for(i=0;i<replicates;i++)
			total[i]=total[i]+node->countTreat[i];
		node=node->Next;
	}
	for(i=0;i<replicates;i++)
		left[i]=ReadNum;

	node=head;
	while(node)
	{
		for(i=0;i<replicates;i++)
		{
			if(total[i]==ReadNum)
				continue;
			node->countTreat[i]=1.0*node->countTreat[i]/total[i]*ReadNum+0.5;
			left[i]=left[i]-node->countTreat[i];
		}
		node=node->Next;
	}

//left to zero
	for(i=0;i<replicates;i++)
	{
		if(total[i]==ReadNum||left[i]==0)
			continue;
		if(left[i]>0)
		{
			head->countTreat[i]=head->countTreat[i]+left[i];
			continue;
		}
		node=head;
		while(node)
		{
			if(node->countTreat[i]+left[i]>=0)
			{
				node->countTreat[i]=node->countTreat[i]+left[i];
				left[i]=0;
				break;
			}
			else
			{
				left[i]=left[i]+node->countTreat[i];
				node->countTreat[i]=0;
				node=node->Next;
			}
		}
	}
		
	if(Flag_control_treat==1)
	{
		for(i=0;i<replicates;i++)
			total[i]=0;
		node=head;
		while(node)
		{
			for(i=0;i<replicates;i++)
				total[i]=total[i]+node->countControl[i];
			node=node->Next;
		}
		for(i=0;i<replicates;i++)
			left[i]=ReadNum;
		node=head;
		while(node)
		{
			for(i=0;i<replicates;i++)
			{
				node->countControl[i]=1.0*node->countControl[i]/total[i]*ReadNum+0.5;
				left[i]=left[i]-node->countControl[i];
			}
			node=node->Next;
		}
	//left to zero
		for(i=0;i<replicates;i++)
		{
			if(total[i]==ReadNum||left[i]==0)
				continue;
			if(left[i]>0)
			{
				head->countControl[i]=head->countControl[i]+left[i];
				continue;
			}
			node=head;
			while(node)
			{
				if(node->countControl[i]+left[i]>=0)
				{
					node->countControl[i]=node->countControl[i]+left[i];
					left[i]=0;
					break;
				}
				else
				{
					left[i]=left[i]+node->countControl[i];
					node->countControl[i]=0;
					node=node->Next;
				}
			}
		}
	}
	free(total);
	free(left);
}

void Abundance_for_Control(struct Genome *head,float more_abundance_per,float less_abundance_per,float maxAbundance,float minAbundance,int readNum)
{
	struct Genome *node;
	int total,more,less,num,i,left;
	float *FC,random;
	
	num=0;
	node=head;
	while(node)
	{
		num++;
		node=node->Next;
	}
	more=1.0*num*more_abundance_per+0.5;
	less=1.0*num*more_abundance_per+0.5;
	FC=(float *)malloc(num*sizeof(float));
	for(i=0;i<num;i++)
		FC[i]=1;
	while(1)
	{
		for(i=0;i<num;i++)
		{
			if(FC[i]!=1)
				continue;
			random=1.0*rand()/RAND_MAX;
			if(random>more_abundance_per+less_abundance_per)
				continue;
			random=1.0*rand()/RAND_MAX;
			if(random<more_abundance_per/(more_abundance_per+less_abundance_per)) //more
			{
				if(more==0)
					continue;
				FC[i]=minAbundance+1.0*(maxAbundance-minAbundance)*rand()/RAND_MAX;
					more--;
			}
			else
			{
				if(less==0)
					continue;
				FC[i]=minAbundance+1.0*(maxAbundance-minAbundance)*rand()/RAND_MAX;
				FC[i]=0-FC[i];
				less--;
			}
			if((more+less)==0)
				break;
		}
		if((more+less)==0)
			break;
	}
	
	i=0;
	total=0;
	node=head;
	while(node)
	{
		if(FC[i]==1)
			node->ReadNum[1]=node->ReadNum[0];
		else if(FC[i]>0)
			node->ReadNum[1]=1.0*node->ReadNum[0]/FC[i]+0.5;
		else
			node->ReadNum[1]=-1.0*node->ReadNum[0]*FC[i]+0.5;
		total=total+node->ReadNum[1];
		i++;
		node=node->Next;
	}

	node=head;
	left=readNum;
	while(node)
	{
		node->ReadNum[1]=1.0*node->ReadNum[1]/total*readNum+0.5;
		left=left-node->ReadNum[1];
		node=node->Next;
	}
	if(left>=0)
		head->ReadNum[1]=head->ReadNum[1]+left;
	else
	{
		node=head;
		while(node)
		{
			if(left==0)
				break;
			if(left+node->ReadNum[1]>=0)
			{
				node->ReadNum[1]=node->ReadNum[1]+left;
				break;
			}
			left=left+node->ReadNum[1];
			node->ReadNum[1]=0;
			node=node->Next;
		}
	}
	free(FC);
	return;
}

void Expression_Control(struct StoreGene *begin,int readNum)
{
	struct StoreGene *node;
	int total,left;

	node=begin;
	while(node)
	{
		total=total+node->countControl[0];
		node=node->Next;
	}
	left=readNum;
	node=begin;
	while(node)
	{
		node->countControl[0]=1.0*node->countControl[0]/total*readNum+0.5;
		left=left-node->countControl[0];
		node=node->Next;
	}
	node=begin;
	while(node)
	{
		if(left==0)
			break;
		if(left>0)
		{
			node->countControl[0]=node->countControl[0]+left;
			break;
		}
		if(node->countControl[0]+left>=0)
		{
			node->countControl[0]=node->countControl[0]+left;
			break;
		}
		else
		{
			left=node->countControl[0]+left;
			node->countControl[0]=0;
			node=node->Next;
		}
	}
	return;
}
int main(int argc,char **argv)
{
	struct Genome *headGenome,*nodeGenome;
	struct GenePos *headGenePos,*nodeGenePos;
	struct InputGene *headInputGene,*beforeInputGene,*nodeInputGene;
	struct DNASequences *headDNASequences,*nodeDNASequences;
	struct StoreGene *headStoremRNA,*headStorerRNA,*nodeStoremRNA,*nodeStorerRNA,*returnStore,*newBeginmRNA,*newBeginrRNA;
	float probability[1001],*FC,express_per,up_expression_per,down_expression_per,maxFC,minFC,pnbinom,*err_rand,type[3],subratio[12],random,qual[50],rRNA;
	float more_abundance_per,less_abundance_per,maxAbundance,minAbundance;
	int i,readNum,flag,*count,j,min,total,num,seed,Flag_pair,Len_read,Len_fragment,sd,Flag_strand_specific,Flag_control_treat,replicates,length_max,Turn_Treat[N],Turn_Control[N],upNum,downNum,cDNA_length;
	int FlagStoremRNA,FlagStorerRNA,FlagNewmRNA,FlagNewrRNA;
	char *cDNA,abundance_file[N],config_file[N],index_file[N],output_prefix[N],gene_file[N],info_file[N],*name;
	FILE *fp,*T1[N],*T2[N],*C1[N],*C2[N];
	time_t time_start,time_end;
	const gsl_rng_type *T;
	gsl_rng *r;

	int get_qv_type();
//default parameters
	readNum=10000;
	Flag_pair=1;
	Len_read=150;
	Len_fragment=350;
	sd=20;
	strcpy(config_file,"simulation.config");
	seed=0;
	Flag_strand_specific=0;
	express_per=0.1;
	Flag_control_treat=0;
	up_expression_per=0.1;
	down_expression_per=0.1;
	maxFC=20;
	minFC=2;
	replicates=-1;
	pnbinom=0.33333;
	min=300;
	rRNA=0.1;
	FlagStoremRNA=0;
	FlagStorerRNA=0;
	more_abundance_per=0.2;
	less_abundance_per=0.2;
	maxAbundance=10;
	minAbundance=2;

	if(argc<6)
	{
		printf("Error! Wrong input.\n");
		usage();
		exit(1);
	}
	time_start=time(NULL);
//get input parameters
	for(i=1;i<argc;)
	{
		if(strcmp(argv[i],"-abundance")==0)
		{
			strcpy(abundance_file,argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-index")==0)
		{
			strcpy(index_file,argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-o")==0)
		{
			strcpy(output_prefix,argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-r")==0)
		{
			readNum=atoi(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-single")==0)
		{
			Flag_pair=0;
			i++;
			continue;
		}
		if(strcmp(argv[i],"-l")==0)
		{
			Len_read=atoi(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-fragment")==0)
		{
			Len_fragment=atoi(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-sd")==0)
		{
			sd=atoi(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-config")==0)
		{
			memset(config_file,'\0',N*sizeof(char));
			strcpy(config_file,argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-seed")==0)
		{
			seed=atoi(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-strandspecific")==0)
		{
			Flag_strand_specific=1;
			i++;
			continue;
		}
		if(strcmp(argv[i],"-express")==0)
		{
			express_per=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-gene")==0)
		{
			strcpy(gene_file,argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-rRNA")==0)
		{
			rRNA=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-diff")==0)
		{
			Flag_control_treat=1;
			i++;
			continue;
		}
		if(strcmp(argv[i],"-up")==0)
		{
			up_expression_per=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-down")==0)
		{
			down_expression_per=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-maxFC")==0)
		{
			maxFC=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-minFC")==0)
		{
			minFC=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-more")==0)
		{
			more_abundance_per=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-less")==0)
		{
			less_abundance_per=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-maxAd")==0)
		{
			maxAbundance=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-minAd")==0)
		{
			minAbundance=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-replicate")==0)
		{
			replicates=atoi(argv[i+1]);
			i=i+2;
			continue;
		}
		if(strcmp(argv[i],"-prob")==0)
		{
			pnbinom=atof(argv[i+1]);
			i=i+2;
			continue;
		}
		else
		{
			printf("Waring: there is no this parameter '%s',please see the usage again!\n",argv[i]);
			i++;
		}
	}

//replicates
	if(replicates==-1&&Flag_control_treat==0)
		replicates=1;
	if(replicates==-1&&Flag_control_treat==1)
		replicates=3;

//check input parameters
	if(check_input(abundance_file,index_file,output_prefix,readNum,Len_fragment,Len_read,config_file,express_per,gene_file,up_expression_per,down_expression_per,minFC,maxFC,pnbinom,rRNA,more_abundance_per,less_abundance_per,maxAbundance,minAbundance)==0)
		exit(1);
//read parameters for sequencing
	fp=fopen(config_file,"r");
	if(fp==NULL)
	{
		printf("Error! Can not open %s file.\n",config_file);
		exit(1);
	}
	err_rand=(float *)malloc(500*50*sizeof(float));
	memset(err_rand,'\0',500*50*sizeof(float));
	length_max=get_qv_type(fp,type,subratio,err_rand);
	if(length_max<Len_read)
	{
		printf("Error! Under th %s configure file, the max read length is %d bp.\n",config_file,length_max);
		printf("  So, can't simulate reads with %d bp.\n",Len_read);
		free(err_rand);
		exit(1);
	}
	
//Expression probability
	for(i=1;i<=1000;i++)
		probability[i]=(float)pow(i,-1.69);
	for(i=2;i<=1000;i++)
		probability[i]=probability[i]+probability[i-1];
	for(i=1;i<=1000;i++)
		probability[i]=probability[i]/probability[1000];
	probability[0]=0;

	for(i=0;i<50;i++)
		qual[i]=pow(10,(i/(-10.0)));

//read files
	headGenome=Abundance(abundance_file,express_per,up_expression_per,down_expression_per,maxFC,minFC,readNum,Flag_control_treat);
	headGenome=readIndex(index_file,headGenome);
	if(headGenome==NULL)
	{
		printf("Error! None organism in abundance table is in index file.\n");
		exit(1);
	}
	ExpressNum(headGenome,express_per,up_expression_per,down_expression_per);
	if(strlen(gene_file)>0)
		ReadGeneList(gene_file,headGenome);
	if(Flag_control_treat==1&&headGenome->ReadNum[1]==-1)
		Abundance_for_Control(headGenome,more_abundance_per,less_abundance_per,maxAbundance,minAbundance,readNum);

//if repeat, initial
	if(replicates>1)
	{
		T=gsl_rng_default;
		r=gsl_rng_alloc(T);
		gsl_rng_set(r,seed);
	//output file
		
		for(i=1;i<=replicates;i++)
		{
			if(Flag_control_treat==1)
			{
				Turn_Treat[i-1]=0;
				Turn_Control[i-1]=0;
				if(Flag_pair==1)
				{
					flag=Create_output(T1,i,output_prefix,"-Treat",1);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
					flag=Create_output(T2,i,output_prefix,"-Treat",2);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
					flag=Create_output(C1,i,output_prefix,"-Control",1);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
					flag=Create_output(C2,i,output_prefix,"-Control",2);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
				}
				else
				{
					flag=Create_output(T1,i,output_prefix,"-Treat",0);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
					flag=Create_output(C1,i,output_prefix,"-Control",0);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
				}
			}
			else
			{
				Turn_Treat[i-1]=0;
				if(Flag_pair==1)
				{
					flag=Create_output(T1,i,output_prefix,NULL,1);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
					flag=Create_output(T2,i,output_prefix,NULL,2);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
				}
				else
				{
					flag=Create_output(T1,i,output_prefix,NULL,0);
					if(flag==0)
					{
						CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
						exit(1);
					}
				}
			}
		}
	}
	else
	{
		if(Flag_control_treat==1)
		{
			Turn_Treat[0]=0;
			Turn_Control[0]=0;
			if(Flag_pair==1)
			{
				flag=Create_output(T1,0,output_prefix,"-Treat",1);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
				flag=Create_output(T2,0,output_prefix,"-Treat",2);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
				flag=Create_output(C1,0,output_prefix,"-Control",1);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
				flag=Create_output(C2,0,output_prefix,"-Control",2);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
			}
			else
			{
				flag=Create_output(T1,0,output_prefix,"-Treat",0);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
				flag=Create_output(C1,0,output_prefix,"-Control",0);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
			}
		}
		else
		{
			Turn_Treat[0]=0;
			if(Flag_pair==1)
			{
				flag=Create_output(T1,0,output_prefix,NULL,1);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
				flag=Create_output(T2,0,output_prefix,NULL,2);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
			}
			else
			{
				flag=Create_output(T1,0,output_prefix,NULL,0);
				if(flag==0)
				{
					CleanMemory(headGenome,err_rand,replicates,NULL,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
					exit(1);
				}
			}
		}
	}
//start rRNA
	if(rRNA>0)  //only those have rRNA fasta file
	{
		nodeGenome=headGenome;
		while(nodeGenome!=NULL)
		{
			if(nodeGenome->MaxLength==0)
			{
				printf("The %s doesn't have sequence information in index file!\n",nodeGenome->Name);
				nodeGenome=nodeGenome->Next;
				continue;
			}
			if(nodeGenome->rRNA==0)
			{
				printf("The %s doesn't have rRNA sequence!\n", nodeGenome->Name);
				nodeGenome=nodeGenome->Next;
				continue;
			}

			FlagNewrRNA=0;
			newBeginrRNA=NULL;
			if(nodeGenome->rRNA_path!=NULL)
			{
				headDNASequences=ReadFastaFile(nodeGenome->rRNA_path,min);
				if(headDNASequences==NULL)
				{
					printf("Error! Can't open %s file of %s genome!\n",nodeGenome->rRNA_path,nodeGenome->Name);
					nodeGenome=nodeGenome->Next;
					continue;
				}
				total=1.0*rRNA*nodeGenome->ReadNum[0]+0.5;
				num=total/nodeGenome->rRNA+0.5;
	                        while(headDNASequences)
				{
					cDNA_length=strlen(headDNASequences->Fasta)+1;
					cDNA=(char *)malloc(cDNA_length*sizeof(char));
					memset(cDNA,'\0',cDNA_length*sizeof(char));
					strcpy(cDNA,headDNASequences->Fasta);

					cDNA_length=strlen(headDNASequences->Name)+1;
					name=(char *)malloc(cDNA_length*sizeof(char));
					memset(name,'\0',cDNA_length*sizeof(char));
					strcpy(name,headDNASequences->Name);
					returnStore=StoreGeneList(Flag_control_treat,1,num,1,cDNA,nodeGenome->Name,name);
					if(FlagStorerRNA==0)
					{
						headStorerRNA=returnStore;
						nodeStorerRNA=headStorerRNA;
						FlagStorerRNA++;
					}
					else
					{
						nodeStorerRNA->Next=returnStore;
						nodeStorerRNA=returnStore;
					}
					if(Flag_control_treat==1&&FlagNewrRNA==0)
					{
						FlagNewrRNA=1;
						newBeginrRNA=returnStore;
					}
					nodeDNASequences=headDNASequences;
					headDNASequences=headDNASequences->Next;
					free(nodeDNASequences->Name);
					free(nodeDNASequences->Fasta);
					free(nodeDNASequences);
					total=total-num;
					if(total<=0)
						break;
				}
				if(total>0&&FlagStorerRNA)
					nodeStorerRNA->countTreat[0]=nodeStorerRNA->countTreat[0]+total;
			}
			if(Flag_control_treat==1&&FlagNewrRNA==1)
				Expression_Control(newBeginrRNA,(1.0*rRNA*nodeGenome->ReadNum[1]+0.5));
			nodeGenome=nodeGenome->Next;
		}
	}
//start mRNA
	nodeGenome=headGenome;
	while(nodeGenome!=NULL)
	{
		if(nodeGenome->MaxLength==0)
		{
			printf("The %s doesn't have sequence information in index file!\n",nodeGenome->Name);
			nodeGenome=nodeGenome->Next;
			continue;
		}
		if(nodeGenome->mRNA==0)
		{
			printf("The %s doesn't have cDNA sequence!\n", nodeGenome->Name);
			nodeGenome=nodeGenome->Next;
			continue;
		}
	//expression
		if(nodeGenome->gene_list!=NULL) //user define the count
		{
			nodeInputGene=nodeGenome->gene_list;
			total=0;
			while(nodeInputGene)
			{
				total+=nodeInputGene->count;
				nodeInputGene=nodeInputGene->Next;
			}
			nodeInputGene=nodeGenome->gene_list;
			num=1.0*(1-rRNA)*nodeGenome->ReadNum[0]+0.5; //delete rRNA
			while(nodeInputGene)
			{
				nodeInputGene->count=(int)(1.0*nodeInputGene->count/total*(1-rRNA)*nodeGenome->ReadNum[0]+0.5);
				num=num-nodeInputGene->count;
				nodeInputGene=nodeInputGene->Next;
			}
			nodeInputGene=nodeGenome->gene_list;
			while(nodeInputGene)
			{
				if(nodeInputGene->count+num>=0)
				{
					nodeInputGene->count=nodeInputGene->count+num;
					break;
				}
				else
				{
					num=nodeInputGene->count+num;
					nodeInputGene->count=0;
					nodeInputGene=nodeInputGene->Next;
				}
			}
		}
		else
		{
			count=(int *)malloc(nodeGenome->mRNA*sizeof(int));
			for(i=0;i<nodeGenome->mRNA;i++)
				count[i]=0;
			i=0;
			total=0;
			while(i<=nodeGenome->DE[0]-1)
			{
				j=(int)(1.0*rand()/RAND_MAX*nodeGenome->mRNA); //from 0 to num_mRNA-1
				if(j>=nodeGenome->mRNA)
					continue;
				if(count[j]>0)
					continue;
				count[j]=GenerateRandomCount(probability);
				total+=count[j];
				i++;
			}
		//normalize
			num=1.0*(1-rRNA)*nodeGenome->ReadNum[0]+0.5;
			for(i=0;i<nodeGenome->mRNA;i++)
			{
				if(count[i]==0)
					continue;
				count[i]=(int)(1.0*count[i]/total*(1-rRNA)*nodeGenome->ReadNum[0]+0.5);
				num=num-count[i];
			}
			i=0;
			while(1)  //adjust the first one
			{
				if(count[i]!=0)
				{
					if(count[i]+num>=0)
					{
						count[i]=count[i]+num;
						break;
					}
					else
					{
						num=num+count[i];
						count[i]=0;
					}
				}
				i++;
			}

		//Fold change
			if(Flag_control_treat==1)
			{
				FC=(float *)malloc(nodeGenome->mRNA*sizeof(float));
				memset(FC,'\0',nodeGenome->mRNA*sizeof(float));
				total=0; //total express genes
				for(i=0;i<nodeGenome->mRNA;i++)
				{
					if(count[i]==0)
						FC[i]=0;
					else
					{
						FC[i]=1;
						total++;
					}
				}
				upNum=(int)nodeGenome->DE[1]; //up-number
				downNum=(int)nodeGenome->DE[2];
				while(1)
				{
					for(j=0;j<nodeGenome->mRNA;j++)
					{
						if(FC[j]==0)
							continue;
						if(FC[j]!=1)
							continue;
						random=1.0*rand()/RAND_MAX;
						if(random>(nodeGenome->DE[1]+nodeGenome->DE[2])/nodeGenome->DE[0])
							continue;
						random=1.0*rand()/RAND_MAX;
						if(random<nodeGenome->DE[1]/(nodeGenome->DE[1]+nodeGenome->DE[2])) //up
						{
							if(upNum==0)
								continue;
							FC[j]=nodeGenome->DE[4]+1.0*(nodeGenome->DE[3]-nodeGenome->DE[4])*rand()/RAND_MAX;
							upNum--;
						}
						else
						{
							if(downNum==0)
								continue;
							FC[j]=nodeGenome->DE[4]+1.0*(nodeGenome->DE[3]-nodeGenome->DE[4])*rand()/RAND_MAX;
							FC[j]=0-FC[j];
							downNum--;
						}
						total--;
						if((upNum+downNum)==0)
							break;
						if(total==0)
							break;
					}
					if((upNum+downNum)==0)
						break;
					if(total==0)
						break;
				}
			}//end-FC
		}
	//cDNA & simulate
		FlagNewmRNA=0;
		FlagNewrRNA=0;
		newBeginmRNA=NULL;
		newBeginrRNA=NULL;
		if(nodeGenome->transcript_path!=NULL)
		{
			headDNASequences=ReadFastaFile(nodeGenome->transcript_path,min);
			if(headDNASequences==NULL)
			{
				printf("Error! Can't open %s file of %s genome!\n",nodeGenome->transcript_path,nodeGenome->Name);
				nodeGenome=nodeGenome->Next;
				continue;
			}
			if(nodeGenome->gene_list!=NULL) //user define
			{
				nodeDNASequences=headDNASequences;
				headInputGene=nodeGenome->gene_list;
				while(nodeDNASequences)
				{
					flag=0;
					if(headInputGene==NULL)
						break;
					if(strcmp(headInputGene->Name,nodeDNASequences->Name)==0)
					{
						nodeInputGene=headInputGene;
						headInputGene=headInputGene->Next;
						flag++;
					}
					else
					{
						beforeInputGene=headInputGene;   //delete the node for search again
						while(beforeInputGene->Next)
						{
							if(strcmp(beforeInputGene->Next->Name,nodeDNASequences->Name)==0)
							{
								flag++;
								nodeInputGene=beforeInputGene->Next;
								beforeInputGene->Next=nodeInputGene->Next;
								break;
							}
							beforeInputGene=beforeInputGene->Next;
						}
					}
					if(flag==0)
					{
						nodeDNASequences=nodeDNASequences->Next;
						continue;
					}
					cDNA_length=strlen(nodeDNASequences->Fasta)+1;
					cDNA=(char *)malloc(cDNA_length*sizeof(char));
					memset(cDNA,'\0',cDNA_length*sizeof(char));
					strcpy(cDNA,nodeDNASequences->Fasta);
					cDNA_length=strlen(nodeInputGene->Name)+1;
					name=(char *)malloc(cDNA_length*sizeof(char));
					memset(name,'\0',cDNA_length*sizeof(char));
					strcpy(name,nodeInputGene->Name);
					if(Flag_control_treat==1)
						returnStore=StoreGeneList(Flag_control_treat,replicates,nodeInputGene->count,nodeInputGene->FC,cDNA,nodeGenome->Name,name);
					else
						returnStore=StoreGeneList(Flag_control_treat,replicates,nodeInputGene->count,1,cDNA,nodeGenome->Name,name);
					if(FlagStoremRNA==0)
					{
						headStoremRNA=returnStore;
						nodeStoremRNA=headStoremRNA;
						FlagStoremRNA++;
					}
					else
					{
						nodeStoremRNA->Next=returnStore;
						nodeStoremRNA=returnStore;
					}
					if(Flag_control_treat==1&&FlagNewmRNA==0)
					{
						FlagNewmRNA=1;
						newBeginmRNA=returnStore;
					}
					nodeDNASequences=nodeDNASequences->Next;
					free(nodeInputGene);
				}
				while(headInputGene)
				{
					printf("Warning: Can't find the sequence of %s gene in %s file!\n",headInputGene->Name,nodeGenome->transcript_path);
					nodeInputGene=headInputGene;
					headInputGene=headInputGene->Next;
					free(nodeInputGene);
				}
				nodeGenome->gene_list=NULL; //clean
			}
			else
			{
				j=0;
				nodeDNASequences=headDNASequences;
				for(i=0;i<nodeGenome->mRNA;i++)
				{
					if(count[i]==0)
						continue;
					flag=0;
					while(nodeDNASequences)
					{
						if(j==i)
						{
							flag=1;
							break;
						}
						nodeDNASequences=nodeDNASequences->Next;
						j++;
					}
					if(flag==0)
					{
						printf("Warning: Can't find the %d-th gene in %s file!\n",i,nodeGenome->transcript_path);
						continue;
					}
					cDNA_length=strlen(nodeDNASequences->Fasta)+1;
					cDNA=(char *)malloc(cDNA_length*sizeof(char));
					memset(cDNA,'\0',cDNA_length*sizeof(char));
					strcpy(cDNA,nodeDNASequences->Fasta);
					cDNA_length=strlen(nodeDNASequences->Name)+1;
					name=(char *)malloc(cDNA_length*sizeof(char));
					memset(name,'\0',cDNA_length*sizeof(char));
					strcpy(name,nodeDNASequences->Name);
					if(Flag_control_treat==1)
						returnStore=StoreGeneList(Flag_control_treat,replicates,count[i],FC[i],cDNA,nodeGenome->Name,name);
					else
						returnStore=StoreGeneList(Flag_control_treat,replicates,count[i],1,cDNA,nodeGenome->Name,name);
					if(FlagStoremRNA==0)
					{
						headStoremRNA=returnStore;
						nodeStoremRNA=headStoremRNA;
						FlagStoremRNA++;
					}
					else
					{
						nodeStoremRNA->Next=returnStore;
						nodeStoremRNA=returnStore;
					}
					if(Flag_control_treat==1&&FlagNewmRNA==0)
					{
						FlagNewmRNA=1;
						newBeginmRNA=returnStore;
					}
				}
				free(count);
				if(Flag_control_treat==1)
					free(FC);
			}
		//free Fasta
			while(headDNASequences)
			{
				nodeDNASequences=headDNASequences;
				headDNASequences=headDNASequences->Next;
				free(nodeDNASequences->Name);
				free(nodeDNASequences->Fasta);
				free(nodeDNASequences);
			}
		}
		else if(nodeGenome->genome_path!=NULL)  //cDNA from annotation file
		{
			headDNASequences=ReadFastaFile(nodeGenome->genome_path,min);
			if(headDNASequences==NULL)
			{
				printf("Error! Can't open %s file!\n",nodeGenome->genome_path);
				nodeGenome=nodeGenome->Next;
				continue;
			}
			if(nodeGenome->AnnotationType==1)
				headGenePos=ReadGenbank(nodeGenome->annotation_path,min);
			else
				headGenePos=ReadGff(nodeGenome->annotation_path,min);

			if(nodeGenome->gene_list!=NULL) //user define
			{
				nodeDNASequences=headDNASequences;
				headInputGene=nodeGenome->gene_list;
				nodeGenePos=headGenePos;
				while(nodeGenePos)
				{
					if(headInputGene==NULL)
						break;
					flag=0;
					if(strcmp(headInputGene->Name,nodeGenePos->Name)==0)
					{
						nodeInputGene=headInputGene;
						headInputGene=headInputGene->Next;
						flag++;
					}
					else
					{
						beforeInputGene=headInputGene;
						while(beforeInputGene->Next)
						{
							if(strcmp(beforeInputGene->Next->Name,nodeGenePos->Name)==0)
							{
								flag++;
								nodeInputGene=beforeInputGene->Next;
								beforeInputGene->Next=nodeInputGene->Next;
								break;
							}
							beforeInputGene=beforeInputGene->Next;
						}
					}
					if(flag==0)
					{
						nodeGenePos=nodeGenePos->Next;
						continue;
					}
				 //genome->cDNA
					flag=0;
					while(nodeDNASequences)
					{
						if(strcmp(nodeDNASequences->Name,nodeGenePos->Genome)==0)
						{
							flag++;
							break;
						}
						nodeDNASequences=nodeDNASequences->Next;
					}
					if(flag==0)
					{
						printf("Warning: Can't find the sequence of %s genome for %s gene!\n",nodeGenePos->Genome,nodeInputGene->Name);
						free(nodeInputGene);
						nodeGenePos=nodeGenePos->Next;
						continue;
					}
					cDNA_length=getLength(nodeGenePos->Pos);
	                                cDNA=(char *)malloc((cDNA_length+1)*sizeof(char));
					memset(cDNA,'\0',(cDNA_length+1)*sizeof(char));
					flag=GeneratecDNA(nodeDNASequences->Fasta,nodeGenePos->Pos,cDNA,cDNA_length);
					if(flag==0)
					{
						printf("Error! Can't generate the cDNA sequence for %s gene!\n",nodeInputGene->Name);
						free(nodeInputGene);
						free(cDNA);
						nodeGenePos=nodeGenePos->Next;
						continue;
					}
					cDNA_length=strlen(nodeInputGene->Name)+1;
					name=(char *)malloc(cDNA_length*sizeof(char));
					memset(name,'\0',cDNA_length*sizeof(char));
					strcpy(name,nodeInputGene->Name);
					if(Flag_control_treat==1)
						returnStore=StoreGeneList(Flag_control_treat,replicates,nodeInputGene->count,nodeInputGene->FC,cDNA,nodeGenome->Name,name);
					else
						returnStore=StoreGeneList(Flag_control_treat,replicates,nodeInputGene->count,1,cDNA,nodeGenome->Name,name);
					if(FlagStoremRNA==0)
					{
						headStoremRNA=returnStore;
						nodeStoremRNA=headStoremRNA;
						FlagStoremRNA++;
					}
					else
					{
						nodeStoremRNA->Next=returnStore;
						nodeStoremRNA=returnStore;
					}
					if(Flag_control_treat==1&&FlagNewmRNA==0)
					{
						FlagNewmRNA=1;
						newBeginmRNA=returnStore;
					}
					free(nodeInputGene);
					nodeGenePos=nodeGenePos->Next;
				}
				while(headInputGene)
				{
					printf("Warning: Can't find the sequence of %s gene in %s file!\n",headInputGene->Name,nodeGenome->annotation_path);
					nodeInputGene=headInputGene;
					headInputGene=headInputGene->Next;
					free(nodeInputGene);
				}
				nodeGenome->gene_list=NULL;//clean
			}
			else
			{
				j=0;
				nodeDNASequences=headDNASequences;
				nodeGenePos=headGenePos;
				for(i=0;i<nodeGenome->mRNA;i++)
				{
					if(count[i]==0)
						continue;
					flag=0;
					while(nodeGenePos)
					{
						if(nodeGenePos->rRNA_flag==1) //not rRNA
						{
							nodeGenePos=nodeGenePos->Next;
							continue;
						}
						if(j==i)
						{
							flag=1;
							break;
						}
						nodeGenePos=nodeGenePos->Next;
						j++;
					}
					if(flag==0)
					{
						printf("Warning: Can't find the %d-th gene in %s file!\n",i,nodeGenome->annotation_path);
						continue;
					}
				//genome->cDNA
					flag=0;
					while(nodeDNASequences)
					{
						if(strcmp(nodeDNASequences->Name,nodeGenePos->Genome)==0)
						{
							flag++;
							break;
						}
						nodeDNASequences=nodeDNASequences->Next;
					}
					if(flag==0)
					{
						printf("Warning: Can't find the sequence of %s genome for %d-th gene!\n",nodeGenePos->Genome,i);
						continue;
					}
					cDNA_length=getLength(nodeGenePos->Pos);
					cDNA=(char *)malloc((cDNA_length+1)*sizeof(char));
					memset(cDNA,'\0',(cDNA_length+1)*sizeof(char));
					flag=GeneratecDNA(nodeDNASequences->Fasta,nodeGenePos->Pos,cDNA,cDNA_length);
					if(flag==0)
					{
						printf("Error! Can't generate the cDNA sequence for %d-th gene!\n",i);
						free(cDNA);
						continue;;
					}
					cDNA_length=strlen(nodeGenePos->Name)+1;
					name=(char *)malloc(cDNA_length*sizeof(char));
					memset(name,'\0',cDNA_length*sizeof(char));
					strcpy(name,nodeGenePos->Name);
					if(Flag_control_treat==1)
						returnStore=StoreGeneList(Flag_control_treat,replicates,count[i],FC[i],cDNA,nodeGenome->Name,name);
					else
						returnStore=StoreGeneList(Flag_control_treat,replicates,count[i],1,cDNA,nodeGenome->Name,name);
					if(FlagStoremRNA==0)
					{
						headStoremRNA=returnStore;
						nodeStoremRNA=headStoremRNA;
						FlagStoremRNA++;
					}
					else
					{
						nodeStoremRNA->Next=returnStore;
						nodeStoremRNA=returnStore;
					}
					if(Flag_control_treat==1&&FlagNewmRNA==0)
					{
						FlagNewmRNA=1;
						newBeginmRNA=returnStore;
					}
				}
				free(count);
				if(Flag_control_treat==1)
					free(FC);
			}
		//rRNA from annotation file
			if(rRNA>0&&nodeGenome->rRNA>0&&(nodeGenome->rRNA_path==NULL))
			{
				total=1.0*rRNA*nodeGenome->ReadNum[0]+0.5;
				num=total/nodeGenome->rRNA+0.5;

				nodeDNASequences=headDNASequences;
				nodeGenePos=headGenePos;
				while(nodeGenePos)
				{
					if(nodeGenePos->rRNA_flag==0)
					{
						nodeGenePos=nodeGenePos->Next;
						continue;
					}
					flag=0;
			//genome->cDNA
					while(nodeDNASequences)
					{
						if(strcmp(nodeDNASequences->Name,nodeGenePos->Genome)==0)
						{
							flag++;
							break;
						}
						nodeDNASequences=nodeDNASequences->Next;
					}
					if(flag==0)
					{
						printf("Warning: Can't find the sequence of %s genome for %s gene!\n",nodeGenome->Name,nodeGenePos->Name);
						nodeGenePos=nodeGenePos->Next;
						continue;
					}
					cDNA_length=getLength(nodeGenePos->Pos);
					cDNA=(char *)malloc((cDNA_length+1)*sizeof(char));
					memset(cDNA,'\0',(cDNA_length+1)*sizeof(char));
					flag=GeneratecDNA(nodeDNASequences->Fasta,nodeGenePos->Pos,cDNA,cDNA_length);
					if(flag==0)
					{
						printf("Error! Can't generate the cDNA sequence for %s gene!\n",nodeGenePos->Name);
						nodeGenePos=nodeGenePos->Next;
						free(cDNA);
						continue;;
					}
					cDNA_length=strlen(nodeGenePos->Name)+1;
					name=(char *)malloc(cDNA_length*sizeof(char));
					memset(name,'\0',cDNA_length*sizeof(char));
					strcpy(name,nodeGenePos->Name);
					returnStore=StoreGeneList(Flag_control_treat,1,num,1,cDNA,nodeGenome->Name,name);
					if(FlagStorerRNA==0)
					{
						headStorerRNA=returnStore;
						nodeStorerRNA=headStorerRNA;
						FlagStorerRNA++;
					}
					else
					{
						nodeStorerRNA->Next=returnStore;
						nodeStorerRNA=returnStore;
					}
					if(Flag_control_treat==1&&FlagNewrRNA==0)
					{
						FlagNewrRNA=1;
						newBeginrRNA=returnStore;
					}
					nodeGenePos=nodeGenePos->Next;
					total=total-num;
					if(total<=0)
						break;
				}
				if(total>0&&FlagStorerRNA)
					nodeStorerRNA->countTreat[0]=nodeStorerRNA->countTreat[0]+total;
			}
		//free Fasta &&GenePos
			while(headDNASequences)
			{
				nodeDNASequences=headDNASequences;
				headDNASequences=headDNASequences->Next;
				free(nodeDNASequences->Name);
				free(nodeDNASequences->Fasta);
				free(nodeDNASequences);
			}
			while(headGenePos)
			{
				nodeGenePos=headGenePos;
				headGenePos=headGenePos->Next;
				free(nodeGenePos->Name);
				free(nodeGenePos->Genome);
				free(nodeGenePos->Pos);
				free(nodeGenePos);
			}
		}
		if(FlagNewmRNA==1)
			Expression_Control(newBeginmRNA,(1.0*(1-rRNA)*nodeGenome->ReadNum[1]+0.5));
		if(FlagNewrRNA==1)
			Expression_Control(newBeginrRNA,(1.0*rRNA*nodeGenome->ReadNum[1]+0.5));
		nodeGenome=nodeGenome->Next;
	}
//end main
	strcpy(info_file,output_prefix);
	strcat(info_file,"-info.txt");
	fp=fopen(info_file,"w");
	if(fp==NULL)
	{
		printf("Error! Can't create the %s file.\n",info_file);
		CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
		exit(1);
	}
	fprintf(fp,"Genome\tGene");
	if(Flag_control_treat==1)
		fprintf(fp,"\tFC");
	if(replicates==1)
		fprintf(fp,"\tTreat");
	else
	{
		for(i=1;i<=replicates;i++)
			fprintf(fp,"\tTreat%d",i);
	}
	if(Flag_control_treat==1)
	{
		if(replicates==1)
			fprintf(fp,"\tControl");
		else
		{
			for(i=1;i<=replicates;i++)
				fprintf(fp,"\tControl%d",i);
		}
	}
	fprintf(fp,"\n");
	//normalize
	if(FlagStorerRNA)
	{
		total=1.0*readNum*rRNA+0.5;
		Normalize(headStorerRNA,1,Flag_control_treat,total);
		while(headStorerRNA)
		{
			fprintf(fp,"%s\trRNA:%s",headStorerRNA->Genome,headStorerRNA->Gene);
			if(Flag_control_treat==1)
			{
				fprintf(fp,"\t1");  //Fold change
				for(i=0;i<replicates;i++)
				{
					fprintf(fp,"\t%d",headStorerRNA->countTreat[0]);
					Turn_Treat[i]=NeSSM(Flag_pair,Flag_strand_specific,headStorerRNA->countTreat[0],Len_read,Len_fragment,sd,headStorerRNA->Fasta,T1[i],T2[i],headStorerRNA->Genome,headStorerRNA->Gene,Turn_Treat[i],type,subratio,err_rand,qual);
				}
				for(i=0;i<replicates;i++)
				{
					fprintf(fp,"\t%d",headStorerRNA->countTreat[0]);
					Turn_Control[i]=NeSSM(Flag_pair,Flag_strand_specific,headStorerRNA->countControl[0],Len_read,Len_fragment,sd,headStorerRNA->Fasta,C1[i],C2[i],headStorerRNA->Genome,headStorerRNA->Gene,Turn_Control[i],type,subratio,err_rand,qual);
				}
			}
			else
			{
				for(i=0;i<replicates;i++)
				{
					fprintf(fp,"\t%d",headStorerRNA->countTreat[0]);
					Turn_Treat[i]=NeSSM(Flag_pair,Flag_strand_specific,headStorerRNA->countTreat[0],Len_read,Len_fragment,sd,headStorerRNA->Fasta,T1[i],T2[i],headStorerRNA->Genome,headStorerRNA->Gene,Turn_Treat[i],type,subratio,err_rand,qual);
				}
			}
			fprintf(fp,"\n");
			nodeStorerRNA=headStorerRNA;
			headStorerRNA=headStorerRNA->Next;
			free(nodeStorerRNA->Gene);
			free(nodeStorerRNA->Fasta);
			free(nodeStorerRNA->countTreat);
			free(nodeStorerRNA);
		}
		readNum=readNum-total;  //mRNA
	}
	if(replicates>1)
		Cal_replicates(Flag_control_treat,replicates,pnbinom,headStoremRNA,r);
	Normalize(headStoremRNA,replicates,Flag_control_treat,readNum);
	while(headStoremRNA)
	{
		fprintf(fp,"%s\t%s",headStoremRNA->Genome,headStoremRNA->Gene);
		if(Flag_control_treat==1)
		{
			fprintf(fp,"\t%f",headStoremRNA->FC);
			for(i=0;i<replicates;i++)
				fprintf(fp,"\t%d",headStoremRNA->countTreat[i]);
			for(i=0;i<replicates;i++)
				fprintf(fp,"\t%d",headStoremRNA->countControl[i]);
		}
		else
		{
			for(i=0;i<replicates;i++)
				fprintf(fp,"\t%d",headStoremRNA->countTreat[i]);
		}
		fprintf(fp,"\n");
		Simulate(Flag_control_treat,replicates,Turn_Treat,Turn_Control,headStoremRNA,Flag_pair,Flag_strand_specific,Len_read,Len_fragment,sd,T1,T2,C1,C2,type,subratio,err_rand,qual);
		nodeStoremRNA=headStoremRNA;
		headStoremRNA=headStoremRNA->Next;
		free(nodeStoremRNA->Gene);
		free(nodeStoremRNA->Fasta);
		free(nodeStoremRNA->countTreat);
		if(Flag_control_treat==1)
			free(nodeStoremRNA->countControl);
		free(nodeStoremRNA);
	}
	CleanMemory(headGenome,err_rand,replicates,r,T1,T2,C1,C2,Flag_control_treat,Flag_pair);
	fclose(fp);

	time_end=time(NULL);
	printf("This simulation takes %f seconds.\n",difftime(time_end,time_start));
	return 0;
}
//======================================================================================
//function:get_qv_type--get error config and return the simulation type
//======================================================================================
int get_qv_type(FILE *fc,float type[],float subratio[],float *err_rand)
{
	int circle_f;
	char qv_config[28000];
	float value[50];
	int length_max=0;
	int get_value();

	while(fgets(qv_config,28000,fc)!=NULL)
	{
		qv_config[strlen(qv_config)-1]='\0';
		if(strstr(qv_config,"illumina_rand")!=NULL)
		{
			get_value(qv_config,value);
			for(circle_f=0;circle_f<50;circle_f++)
				err_rand[length_max*50+circle_f]=value[circle_f];
			length_max++;
			continue;
		}
		if(strstr(qv_config,"illumina_type")!=NULL)
		{
			get_value(qv_config,type);
			continue;
		}

		if(strstr(qv_config,"illumina_sub")!=NULL)
		{
			get_value(qv_config,subratio);
			continue;
		}
	}
	fclose(fc);
	return length_max;
}
///======================================================================================
//function:get_value
//======================================================================================
int get_value(char qv_config[],float need[])
{
	int circle_f,length;
	char *part[2],*value[4000];
	split(part,qv_config,"=");
	length=split(value,part[1],":");
	for(circle_f=0;circle_f<length;circle_f++)
	{
		need[circle_f]=atof(value[circle_f]);
	}
	return length;
}
//======================================================================================
//function:quality--return an error value
//======================================================================================
int quality(int pos_base,float *err_rand)  //pos_base:from zero
{
	float random;
	int circle,step,turn;
	step=50;
	random=(float)(rand())/(float)(RAND_MAX);

	if(random<=err_rand[step*pos_base+30])
	{
		for(circle=(step*pos_base+30);circle>step*pos_base;circle--)
		{
			if((random>err_rand[circle-1])&&(random<=err_rand[circle]))
			{
				turn=circle-step*pos_base;
				return turn;
			}
		}
		return 0;
	}
	else
	{
		for(circle=(step*pos_base+31);circle<step*(pos_base+1);circle++)
		{
			if((random>err_rand[circle-1])&&(random<=err_rand[circle]))
			{
				turn=circle-step*pos_base;
				return turn;
			}
		}
	}
	return 49;
}
