//	extract the keyword=value (s) from an input string 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define BIG_BUF	255
#define SEPS " \t"
#define MAXKEY	20

int main( int argc, char *argv[])
{
	int quietmode=0;
	int tabs=0;
	int len;
	char buf[BIG_BUF],*key[MAXKEY],*newline,value[20];
	char *def=0;
	int nkey=0;

	if(argc < 2 ){
		fprintf(stderr,"usage: getkey [quiet] [tabs] [default=] keyword ...\n");
		exit(1);
	}
	
	for(int i=1; i<argc; i++){
		if(strcmp(argv[i],"quiet")== 0){
			quietmode=1;
			continue;
		}
		if(strcmp(argv[i],"tabs")== 0){
			tabs=1;
			continue;
		}
		if(strncmp(argv[i],"default=",8)==0){
			quietmode=1;
			def=argv[i]+8;
			continue;
		}
		sprintf(buf,"%s=",argv[i]);
		key[nkey++]=strdup(buf);;
		if(nkey==MAXKEY){
			fprintf(stderr,"getkey:  too many keys\n");
			exit(1);
		}
	}

	if(nkey==0){
		fprintf(stderr,"getkey: no keyword specified\n");
		exit(1);
	}

	while(fgets(buf,BIG_BUF,stdin)){
		newline=strstr(buf,"\n");
		if(newline)
			*newline=0;
		for(int i=0; i<nkey; i++){
			len=strlen(key[i]);
			char *gotit=strstr(buf,key[i]);
			if(gotit==NULL){
				if(def)
					printf("%s\n",def);
				if(!quietmode)
					fprintf(stderr,"getkey: keyword %s not found in '%s'\n",
							key[i],buf);
				exit(1);
			}else{
				sscanf(gotit+len,"%s",value);
				printf("%s%s",value,tabs?"\t":" ");
			}
		}
		printf("\n");
	}
}
