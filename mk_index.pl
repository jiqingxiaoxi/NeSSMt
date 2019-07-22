use strict; use warnings;
use Getopt::Long;

my $check;
my $genome_list;
my $annotation_list;
my $transcript_list;
my $rRNA_list;
my $flag;
my $i;
my @list;
my @store;
my $turn=0;
my %hash;
my $output;
my $accession;
my $pos;
my $name;
my $min=300;

my $help="Usage:
  perl $0 --transcript|-t <transcript>  -output|-o <output_file> [options]*

Arguments:
  --transcript|-t <string>
    the transcript files, or the directories containing the transcripts files.
    files or directories are separated by comma.
    if input the directories of transcript files, neither genome nor rRNA gene files are allowed in those directories.
    transcript files must be suffixed with ".fa", ".fasta" or ".fas".
  --genome|-g <string>
    the genome files, or the directories containing the genome files.
    files or directories are separated by comma.
    if input the directories of genome files, neither transcript nor rRNA gene files are allowed in those directories.
    genome files must be suffixed with ".fa", ".fasta" or ".fas".
    it must be used with --annotation|-a argument(below).
  --annotation|-a <string>
    the annotation files, or the directories containing the annotation files.
    files or directories are separated by comma.
    annotation files must be suffixed with ".genbank", ".gb", ".gbff" or ".gff3".
    it must be used with --genome|-g argument(above).
  --rRNA|-r <string>
    the rRNA gene files, or the directories containing the rRNA gene files.
    files or directories are separated by comma.
    if input the directories of rRNA gene files, neither transcript nor genome files are allowed in those directories.
    rRNA gene files must be suffixed with ".fa", ".fasta" or ".fas".
  --output|-o <string>
    the output file.
  --min <int>
    the min length of transcript or rRNA gene sequence used to simulate reads.
    default: 300.
  --help|-h
    print help information.\n";

$check=GetOptions("genome|g=s"=>\$genome_list,"annotation|a=s"=>\$annotation_list,"transcript|t=s"=>\$transcript_list,"rRNA|r=s"=>\$rRNA_list,"output|o=s"=>\$output,"help|h"=>\$flag,"min=i"=>\$min);
if($check==0||$flag)
{
	print $help;
	exit;
}
if($genome_list&&(!$annotation_list))
{
	print "Error! The --genome|-g argument must be with --annotation|-a!\n\n";
	print $help;
	exit;
}
if($annotation_list&&(!$genome_list))
{
	print "Error! The --annotation|-a argument must be with --genome|-a!\n\n";
	print $help;
	exit;
}

if(!$output)
{
	print "Error! Don't have the --output argument!\n\n";
	print $help;
	exit;
}

if($transcript_list)
{
	@list=split(",",$transcript_list);
	for($i=0;$i<@list;$i++)
	{
		if(-f $list[$i])
		{
			$accession=getAccession($list[$i]);
			$hash{$accession}=$turn;
			setdefault($turn);
			$store[$turn][0]=$accession;
			$store[$turn][5]=$list[$i];
			$store[$turn][3]=readFasta($list[$i],$turn,$min);
			$turn++;
			next;
		}
		if(-d $list[$i])
		{
			$turn=scanFasta($list[$i],$turn,1,$min);
		}
	}
}

if($rRNA_list)
{
	@list=split(",",$rRNA_list);
	for($i=0;$i<@list;$i++)
	{
		if(-f $list[$i])
                {
			$accession=getAccession($list[$i]);
			if(exists $hash{$accession})
			{
				$pos=$hash{$accession};
				$store[$pos][9]=$list[$i];
				$store[$pos][4]=readFasta($list[$i],$pos,$min);
			}
			else
			{
				$hash{$accession}=$turn;
				setdefault($turn);
				$store[$turn][0]=$accession;
				$store[$turn][9]=$list[$i];
				$store[$turn][4]=readFasta($list[$i],$turn,$min);
				$turn++;
				next;
			}
		}
		if(-d $list[$i])
		{
			$turn=scanFasta($list[$i],$turn,2,$min);
		}
	}
}
if($genome_list)
{
	@list=split(",",$genome_list);
	for($i=0;$i<@list;$i++)
	{
		if(-f $list[$i])
		{
			$accession=getAccession($list[$i]);
			if(exists $hash{$accession})
			{
				$pos=$hash{$accession};
				readFasta($list[$i],$pos,$min);
				open(IN,"<$list[$i]") or die "Can't open $list[$i] file!\n";
                                $name=<IN>;
                                close IN;
                                chomp $name;
                                ($name)=$name=~/\s+(.+)$/;
                                while($name=~/\,/)
                                {
                                        ($name)=$name=~/^(.+)\,/;
                                }
                                $store[$pos][1]=$name;
                                $store[$pos][6]=$list[$i];
                        }
                        else
                        {
                                $hash{$accession}=$turn;
                                setdefault($turn);
				readFasta($list[$i],$turn,$min);
                                $store[$turn][0]=$accession;
                                $store[$turn][6]=$list[$i];
				open(IN,"<$list[$i]") or die "Can't open $list[$i] file!\n";
				$name=<IN>;
				close IN;
				chomp $name;
				($name)=$name=~/\s+(.+)$/;
				while($name=~/\,/)
				{
					($name)=$name=~/^(.+)\,/;
				}
				$store[$turn][1]=$name;
                                $turn++;
                                next;
                        }
                }
                if(-d $list[$i])
                {
                        $turn=scanFasta($list[$i],$turn,0,$min);
                }
        }
}
if($annotation_list)
{
	@list=split(",",$annotation_list);
        for($i=0;$i<@list;$i++)
        {
                if(-f $list[$i])
                {
                        $accession=getAccession($list[$i]);
                        if(exists $hash{$accession})
                        {
                                $pos=$hash{$accession};
				if($store[$pos][3]>0&&$store[$pos][4]>0) ##have transcript ref and rRNA ref
				{
					next;
				}
				if($list[$i]=~/gb$/i||$list[$i]=~/genbank$/i||$list[$i]=~/gbff$/i||$list[$i]=~/gbk$/i)
				{
					readGenbank($list[$i],$pos,$min);
				}
				elsif($list[$i]=~/gff3$/i||$list[$i]=~/gff$/i)
				{
					readGff3($list[$i],$pos,$min);
				}
				else
				{
					print "Program can identify the file formate of $list[$i]!\n";
				}
                        }
                        else
                        {
                                $hash{$accession}=$turn;
                                setdefault($turn);
				if($list[$i]=~/gb$/i||$list[$i]=~/genbank$/i||$list[$i]=~/gbff$/i||$list[$i]=~/gbk$/i)   
                                {
                                        readGenbank($list[$i],$turn,$min);
                                }
                                elsif($list[$i]=~/gff3$/i||$list[$i]=~/gff$/i)
                                {
                                        readGff3($list[$i],$turn,$min);
                                }
                                else
                                {
                                        print "Program can identify the file formate of $list[$i]!\n";
					next;
                                }
				$store[$turn][0]=$accession;
                                $turn++;
                                next;
                        }
                }
                if(-d $list[$i])
                {
                        $turn=scanAnnotation($list[$i],$turn,$min);
                }
        }
}

open(OUT,">$output") or die "Can't create $output\n";
for($i=0;$i<$turn;$i++)
{
	print OUT $store[$i][0];
	if($store[$i][1] ne "")
	{
		print OUT ":$store[$i][1]";
	}
	print OUT "\t$store[$i][2]\t$store[$i][3]\t$store[$i][4]";
	if($store[$i][5] ne "")
	{
		print OUT "\tT:$store[$i][5]";
	}
	if($store[$i][9] ne "")
	{
		print OUT "\tR:$store[$i][9]";
	}
	if($store[$i][6] ne "")
	{
		print OUT "\tG:$store[$i][6]";
	}
	if($store[$i][7] ne "")
	{
		print OUT "\t$store[$i][8]\:$store[$i][7]";
	}
	print OUT "\n";
}
close OUT;

sub getAccession
{
	my ($in)=@_;
	while($in=~/\//)
	{
		($in)=$in=~/\/(.+)$/;
	}
	($in)=$in=~/^(.+)\.[\w\d]+$/;
	return $in;
}
sub setdefault
{
	my ($in)=@_;
	$store[$in][0]="";
	$store[$in][1]="";
	$store[$in][2]=0; ##max length
	$store[$in][3]=0; ##num of transcripts
	$store[$in][4]=0; ##num of rRNA
	$store[$in][5]=""; ##transcript_path
	$store[$in][6]=""; ##genome_path
	$store[$in][7]=""; ##annotation_path
	$store[$in][8]=""; ##file formate of annotation
	$store[$in][9]=""; ##rRNA path
}

sub readFasta
{
	my ($file,$turn,$min)=@_;
	my $line;
	my $length=-1;
	my $max=0;
	my $num=0;

	open(IN,"<$file") or die "Can't open $file file!\n";
	while(<IN>)
	{
		chomp;
		$line=$_;
		if($line=~/^>/)
		{
			if($length>=$min)
			{
				$num++;
			}
			if($length>$max)
			{
				$max=$length;
			}
			$length=0;
			next;
		}
		$length=$length+length($line);
	}
	close IN;
	if($length>=$min)
	{
		$num++;
	}
	if($length>$max)
	{
		$max=$length;
	}
	if($max>$store[$turn][2])
	{
		$store[$turn][2]=$max;
	}
	return $num;
}

sub scanFasta
{
	my ($source,$turn,$flag,$min)=@_; ##$flag=0:genome,1:transcript,2:rRNA
    	opendir(DIR,$source);
    	my @dir = readdir(DIR);
	my $i;
	my $num;
	my $accession;
	my $pos;
	my $name;
	my $length;

    	if($source!~/\/$/)
	{
		$source.="/";
    	}
    	for(my $i=0;$i<@dir;$i++)
	{
		if($dir[$i]=~/^\./)
		{
	    		next;
		}
		$dir[$i]=$source.$dir[$i];
		if(-d $dir[$i])
		{
	    
	    		$turn=scanFasta($dir[$i],$turn,$flag,$min);
		}
		elsif(-f $dir[$i])
		{
	    		if($dir[$i]=~/fna$/i||$dir[$i]=~/fa$/i||$dir[$i]=~/fasta$/i||$dir[$i]=~/fas$/i)
			{
				$accession=getAccession($dir[$i]);
				if(exists $hash{$accession})
				{
					$pos=$hash{$accession};
					$num=readFasta($dir[$i],$pos,$min);
					if($flag==0)
					{
						$store[$pos][6]=$dir[$i];
						open(IN,"<$dir[$i]") or die "Can't open $dir[$i] file!\n";
						$name=<IN>;
						close IN;
						chomp $name;
						($name)=$name=~/\s+(.+)$/;
						while($name=~/\,/)
						{
							($name)=$name=~/^(.+)\,/;
						}
						$store[$pos][1]=$name;
					}
					elsif($flag==1)
					{
						$store[$pos][3]=$num;
						$store[$pos][5]=$dir[$i];
					}
					else
					{
						$store[$pos][4]=$num;
						$store[$pos][9]=$dir[$i];
					}
				}
				else
				{
					$hash{$accession}=$turn;
					setdefault($turn);
					$store[$turn][0]=$accession;
					$num=readFasta($dir[$i],$turn,$min);
					if($flag==0)
					{
						$store[$turn][6]=$dir[$i];
						open(IN,"<$dir[$i]") or die "Can't open $dir[$i] file!\n";
						$name=<IN>;
						close IN;
						chomp $name;
						($name)=$name=~/\s+(.+)$/;
						while($name=~/\,/)
						{
							($name)=$name=~/^(.+)\,/;
						}
						$store[$turn][1]=$name;
					}
					elsif($flag==1)
                                        {
						$store[$turn][3]=$num;
						$store[$turn][5]=$dir[$i];
					}
					else
					{
						$store[$turn][4]=$num;
						$store[$turn][9]=$dir[$i];
					}
					$turn++;
				}	
	    		}
		}
		else
		{
	    		next;
		}
    	}
    	close DIR;
	return $turn;
}

sub getLength
{
	my ($in)=@_;
	my @array;
	my $i;
	my $start;
	my $end;
	my $length=0;

	@array=split(",",$in);
	for($i=0;$i<@array;$i++)
	{
		($start,$end)=$array[$i]=~/^(.+)\.\.(.+)$/;
		while($start!~/^\d/)
		{
			($start)=$start=~/(\d+.*)$/;
		}
		while($start!~/\d$/)
		{
			($start)=$start=~/^(\d+)/;
		}
		while($end!~/^\d/)         
                {        
                        ($end)=$end=~/(\d+.*)$/;         
                }                
                while($end!~/\d$/)                        
                {
                        ($end)=$end=~/^(\d+)/;          
                }
		$length=$length+$end-$start+1;
	}
	return $length;
}	

sub readGenbank
{
	my ($file,$turn,$min)=@_;
	my $mRNA=0;
	my $rRNA=0;
	my $line;
	my $flag=0;
	my $id;
	my $yes=0;
	my @gene;
	my $num=0;
	my $i;
	my $length=0;
	my $RNA_flag=0;
	my $CDS_flag=0;
	my $pos_flag;
	my $pos;
	my $PosString;

	open(IN,"<$file") or die "Can't open $file file!\n";
	while(<IN>)
	{
		chomp;
		$line=$_;
		if($line=~/^\s+gene\s\s\s+/&&$line=~/\d+\.\./)
		{
			$flag=1;
			next;
		}
		if($line=~/^\s+\/pseudo$/&&$flag)
		{
			$flag=0;
			$num--;
			next;
		}
	##gene name
		if($line=~/^\s+\/gene=/&&$flag==1)
		{
			($id)=$line=~/\=\"(.+)\"/;
			$gene[$num]=$id;
			$flag++;
			$num++;
			next;
		}
		if($line=~/^\s+\/locus_tag=/&&$flag==1)
		{
			($id)=$line=~/\=\"(.+)\"/;
                        $gene[$num]=$id;
                        $flag++;
			$num++;
                        next;
                }
	##gene is mRNA, CDS or ncRNA
		if(($line=~/^\s+mRNA\s\s\s+/||$line=~/^\s+ncRNA\s\s\s+/||$line=~/^\s+CDS\s\s\s+/)&&$line=~/\d+\.\./)
		{
			$yes=1;
			($PosString)=$line=~/^\s+\w+\s+([\w\d\.\>\<\)\(\,]+)$/;
			$pos_flag=1;
			next;
		}
	##pos string
		if($pos_flag&&($line!~/^\s+\//))
		{
			($pos)=$line=~/^\s+([\w\d\.\>\<\)\(\,]+)$/;
			$PosString=$PosString.$pos;
			next;
		}
		if($pos_flag&&($line=~/^\s+\//))
		{
			$pos_flag=0;
		}
		if($line=~/^\s+\/gene=/&&$yes)
                {
                        ($id)=$line=~/\=\"(.+)\"/;
			$yes=0;
			if($num==0)  ##wrong mRNA
			{
				next;
			}
			if($gene[$num-1] eq $id)
			{
				$length=getLength($PosString);
				if($length>=$min)
				{
					$mRNA++;
				}
				$num--;
			}
			else
			{
				for($i=0;$i<$num;$i++)
				{
					if($gene[$i] eq $id)
					{
						$length=getLength($PosString);
						if($length>=$min)
						{
							$mRNA++;
						}
						$gene[$i]="NULL";
						last;
					}
				}
				##if i!=$num: wrong mRNA
			}
                        next;
                }
                if($line=~/^\s+\/locus_tag=/&&$yes)
                {
			($id)=$line=~/\=\"(.+)\"/;
                        $yes=0;
                        if($num==0)  ##wrong mRNA 
                        {
                                next;
                        }
                        if($gene[$num-1] eq $id)
                        {
				$length=getLength($PosString);
				if($length>=$min)
				{
					$mRNA++;
				}
                                $num--;
                        }
                        else
                        {
                                for($i=0;$i<$num;$i++)
                                {
                                        if($gene[$i] eq $id)
                                        {
						$length=getLength($PosString);
						if($length>=$min)
						{
							$mRNA++;
						}
                                                $gene[$i]="NULL";
                                                last;
                                        }
                                }               
                                ##if i!=$num: wrong mRNA
                        }               
                        next;
                }
	##rRNA
		if($line=~/^\s+rRNA\s\s\s+/&&$line=~/\d+\.\.\d+/)
		{
			$flag=0;
			($PosString)=$line=~/^\s+\w+\s+([\w\d\.\>\<\)\(\,]+)$/;
			$length=getLength($PosString);
			if($length>=$min)
			{
				$rRNA++;
			}
			next;
		}
	##ensure the right $id
		if($line=~/^\s+[\w\-]+\s\s\s+/&&$line=~/\d+\.\./)
		{
			$flag=0;
			$yes=0;
		}
	}
	close IN;
	if($store[$turn][3]==0)
	{
		$store[$turn][3]=$mRNA;
	}
	if($store[$turn][4]==0)
	{
		$store[$turn][4]=$rRNA;
	}
	$store[$turn][7]=$file;
	$store[$turn][8]="GB";
}


sub readGff3
{
	my ($file,$turn,$min)=@_;
	my $line;
	my @array;
	my $mRNA=0;
	my $rRNA=0;
	my $flag=0;
	my $id;
	my $length=0;
	my $CDS_flag=0;
	my $RNA_flag=0;

	open(IN,"<$file") or die "Can't open $file file!\n";
	while(<IN>)
	{
		$line=$_;
		if($line=~/^\#/)
		{
			next;
		}
		@array=split("\t",$line);
	##some line may be empty
		if(@array<9)
		{
			next;
		}
		if($array[2] eq "gene")
                {
			$CDS_flag=0;
			$RNA_flag=0;
			if($length>=$min)
			{
				$mRNA++;
			}
		##gene name
			($id)=$line=~/Name\=(.+)$/;
			while($id=~/\;/)
			{
				($id)=$id=~/^(.+)\;/;
			}
			$length=0;
                        $flag=1;
                        next;
                }
		if($array[2] eq "pseudogene")
		{
			$flag=2;
			$CDS_flag=0;
			$RNA_flag=0;
			next;
		}
                if($flag==1&&($line=~/gbkey=mRNA/||$line=~/gbkey=ncRNA/||$line=~/gbkey=CDS/))
                {
                        $flag++;
			$length=0;
			if($line=~/gbkey=CDS/)
			{
				$CDS_flag=1;
				$length=$array[4]-$array[3]+1;
			}
			else
			{
				$RNA_flag=1;
			}
                        next;
                }
		if($array[2] eq "exon")
		{
			if($RNA_flag)
			{
				$length=$length+$array[4]-$array[3]+1;
				next;
			}
			next;
		}
		if($array[2] eq "CDS")
		{
			if($CDS_flag)
			{
				$length=$length+$array[4]-$array[3]+1;
				next;
			}
			next;
		}

	##rRNA
                if($array[2] eq "rRNA")
                {
			if($length>$min)
                        {
                                $mRNA++;
                        }
                        $flag=2;
			$length=$array[4]-$array[3]+1;
			if($length>=$min)
			{
				$rRNA++;
			}
			$length=0;
			$RNA_flag=0;
			$CDS_flag=0;
                        next;
                }
	##reset the flag
		$RNA_flag=0;
		$CDS_flag=0;
        }
        close IN;
	if($length>=$min)
	{
		$mRNA++;
	}
        if($store[$turn][3]==0)
        {
                $store[$turn][3]=$mRNA;
        }
        if($store[$turn][4]==0)
        {
                $store[$turn][4]=$rRNA;
        }
        $store[$turn][7]=$file;
        $store[$turn][8]="GF";
}

sub scanAnnotation
{
        my ($source,$turn,$min)=@_;
        opendir(DIR,$source);
        my @dir = readdir(DIR);
        my $i;
        my $accession;
        my $pos;
        my $length;

        if($source!~/\/$/)
        {
                $source.="/";
        }
        for(my $i=0;$i<@dir;$i++)
        {
                if($dir[$i]=~/^\./)
                {
                        next;
                }
                $dir[$i]=$source.$dir[$i];
                if(-d $dir[$i])
                {
                        $turn=scanAnnotation($dir[$i],$turn,$min);
                }
                elsif(-f $dir[$i])
                {
                        if($dir[$i]=~/gb$/i||$dir[$i]=~/genbank$/i||$dir[$i]=~/gbff$/i||$dir[$i]=~/gbk$/i)
                        {
                                $accession=getAccession($dir[$i]);
                                if(exists $hash{$accession})
                                {
                                        $pos=$hash{$accession};
                                        readGenbank($dir[$i],$pos,$min);
                                }
                                else
                                {
                                        $hash{$accession}=$turn;
                                        setdefault($turn);
					$store[$turn][0]=$accession;
                                        readGenbank($dir[$i],$turn,$min);
                                        $turn++;
                                }
				next;
                        }
			if($dir[$i]=~/gff3$/i||$dir[$i]=~/gff$/i)
			{
				$accession=getAccession($dir[$i]);
                                if(exists $hash{$accession})
                                {
                                        $pos=$hash{$accession};
                                        readGff3($dir[$i],$pos,$min);   
                                }
                                else
                                {
                                        $hash{$accession}=$turn;
                                        setdefault($turn);
					$store[$turn][0]=$accession;
                                        readGff3($dir[$i],$turn,$min);   
                                        $turn++;
                                }
                                next;
                        }
                }
                else
                {
                        next;
                }
        }
        close DIR;
        return $turn;
}
