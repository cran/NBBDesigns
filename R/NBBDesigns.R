#' NBB Design  of Type 1 (NBBD 1)
#'@name nbbd1
#' @param v  Number of treatments (a prime number)
#'@description When v(>3) is a prime number (prime number is a number which is divisible by itself and by 1) then it will generate a class of neighbour balanced block designs. This class of design is generated using the method of constructing Type 3 designs of (Azais, 1993). It also gives the parameters of the design, information matrix for estimating the contrast pertaining to direct and neighbour effects (both left and right) of the treatments.
#' @return It generates neighbour balanced blocks designs for a given number of treatments (v>3), when v is prime number.
#' @export
#'@references
#'(i) Azais, J.M. (1987)<DOI: 10.1111/j.2517-6161.1987.tb01704.x>."Design of experiments for studying intergenotypic competition";
#'(ii) Azais, J.M., Bailey, R.A. and Monod, H. (1993)<DOI: 10.2307/2532269>."A catalogue of efficient neighbour designs with border plots";
#'(iii) Rees, D. H. (1967)<DOI: 10.2307/2528428>."Some designs of use in serology";
#'
#'
#' @examples
#' library(NBBDesigns)
#' nbbd1(5)
#'
nbbd1<-function(v){
  v= #Number of treatments
    np<-v

  type=1 #Type of Design

  i=2
  while(i<=(np/2)){
    if(np%%i!=0){
      i=i+1
    }else{
      #print(c( "Entered number is not a prime number "), quote=FALSE)

      break
    }
  }
  if(i>np/2 && v>3 ){

    ########################

    mat<-t(matrix(0:(v-1),nrow=v,ncol=v))
    i=0
    while(i<v){
      mat[i+1,]<-(mat[i+1,]*(i+1))
      x1<-mat
      i=i+1

    }
    x1<-(x1%%v)+1
    x1<-x1[-v,]
    matt<-x1
    c1<-matrix("|",nrow=v-1,ncol=1)
    x1<-cbind(x1[,v],c1,x1,c1,x1[,1])

    print("NBBD with left and right border plots", quote=FALSE)
    prmatrix(x1,rowlab=rep("",v-1),collab=rep("",v+4), quote = FALSE)
    #number_of_blocks
    b<-(v-1)

    print(c("Number of blocks",b), quote = FALSE)
    #number_of_replication
    r<-v-1
    print(c("Number of replications",r), quote = FALSE)
    #number_of_treatments
    print(c("Number of treatments",v), quote = FALSE)
    # block_size
    k<-v
    print(c("Block size",k), quote = FALSE)
    #both_side_neighbour_appears
    u<-1
    print(c(" Number of times each treatment appearing as both left and right neighbour to each other treatments",u) ,quote= FALSE)

    ##################


    mat1<-matt
    mat2<-cbind(mat1[,ncol(mat1)],mat1[,c(1:(ncol(mat1)-1))])
    mat3<-cbind(mat1[,c(2:(ncol(mat1)),1)])
    ################################
    incident1<-matrix(0,nrow=length(mat1), ncol=v)
    i=1
    while(i<=v){
      x2<-c(which(t(mat1)==i))
      for(j in x2){
        incident1[j,i]<-(incident1[j,i]+1)
      }
      i=i+1
    }

    ################################
    incident2<-matrix(0,nrow=length(mat1), ncol=v)
    i=1
    while(i<=v){
      x2<-c(which(t(mat2)==i))
      for(j in x2){
        incident2[j,i]<-incident2[j,i]+1
      }
      i=i+1
    }

    ##############################
    incident3<-matrix(0,nrow=length(mat1), ncol=v)
    i=1
    while(i<=v){
      x2<-c(which(t(mat3)==i))
      for(j in x2){
        incident3[j,i]<-incident3[j,i]+1
      }
      i=i+1
    }

    #############################
    #D_matrix

    k=1
    d_mat<-matrix(,nrow=0,ncol=b)
    while(k<=b){
      xd<-matrix(,nrow=length(mat1[k,]),ncol=0)
      id<-matrix(0,nrow=length(mat1[k,]),ncol=b)
      id[,k]=1
      xd<-id
      d_mat<-rbind(d_mat,(xd))


      k=k+1

    }


    ##################################
    x1_mat<-cbind(incident1,incident2,incident3)
    vec1n<-matrix(1,nrow=nrow(incident3),ncol=1)
    x2_mat<-cbind(vec1n,d_mat)
    #################x1 prime x1
    x1_prime_x1<-t(x1_mat)%*% x1_mat
    #################x1 prime x2
    x1_prime_x2<-t(x1_mat)%*% x2_mat
    #################x2 prime x2
    x2_prime_x2<-t(x2_mat)%*% x2_mat

    #########################





    #joint C matrix
    Cmatrix<-(x1_prime_x1)-(x1_prime_x2)%*%MASS::ginv(x2_prime_x2)%*%(t(x1_prime_x2))
    Cmatrix<-round(Cmatrix,digits = 3)
    #print("Joint C matrix", quote=FALSE)
    #print(Cmatrix)


    #The information matrix for estimating the direct effects of treatments
    c11<-Cmatrix[c(1:v),c(1:v)]
    c12<-Cmatrix[c(1:v),c((v+1):(3*v))]
    c21<-t(c12)
    c22<-Cmatrix[c((v+1):(3*v)),c((v+1):(3*v))]

    C_tau<-c11-c12%*%MASS::ginv(c22)%*%c21
    print("The information matrix for estimating the contrast pertaining to the direct effects of treatments", quote=FALSE)
    print(round(C_tau,digits = 3))
    #The information matrix for estimating the left neighbour effects of treatments

    c11<-Cmatrix[c((v+1):(2*v)),c((v+1):(2*v))]
    c121<-Cmatrix[c((v+1):(2*v)),c((1:v))]
    c122<-Cmatrix[c((v+1):(2*v)),c((2*v+1):(3*v))]
    c12<-cbind(c121,c122)
    c21<-t(c12)
    c221<-Cmatrix[c(1:v),c(1:v)]
    c222<-Cmatrix[c(1:v),c(((2*v)+1):(3*v))]
    c223<-Cmatrix[c(((2*v)+1):(3*v)),c(1:v)]
    c224<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
    c225<-cbind(c221,c222)
    c226<-cbind(c223,c224)
    c22<-rbind(c225,c226)
    C_ro<-c11-c12%*%MASS::ginv(c22)%*%c21
    print("The information matrix for estimating the contrast pertaining to the left neighbour effects of treatments", quote=FALSE)
    print(round(C_ro,digits = 3))
    #The information matrix for estimating the right neighbour effects of treatments
    c11<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
    c12<-Cmatrix[c(((2*v)+1):(3*v)),c(1:(2*v))]
    c21<-t(c12)
    c22<-Cmatrix[c(1:(2*v)),c(1:(2*v))]

    C_del<-c11-c12%*%MASS::ginv(c22)%*%c21
    print("The information matrix for estimating the contrast pertaining to the right neighbour effects of treatments", quote=FALSE)
    print(round(C_del,digits = 3))

  }else{
    print("Please enter a correct value",quote=FALSE)
  }
}


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



#' NBB Design  of Type 2 (NBBD 2)
#'@name nbbd2
#' @param v Number of treatments (v is odd number but not a power of a prime number)
#'@description When v is odd number but not a power of a prime number (v<= 100) , the function generates a class of neighbour balanced block designs for a given number of treatments. This design listed as type 1 design (Azais, 1993) which is constructed using the method of Rees (1967) and Azais (1987). It also gives the parameters of the design, information matrix for estimating the contrast pertaining to direct and neighbour effects (both left and right) of the treatments.
#' @return It generates neighbour balanced block designs for a given number of treatments (v>3),  v is odd number but not a power of a prime number.
#' @export
#'@references Azais, J.M., Bailey, R.A. and Monod, H. (1993)<DOI: 10.2307/2532269>."A catalogue of efficient neighbour designs with border plots".
#' @examples
#' library(NBBDesigns)
#' nbbd2(9)
nbbd2<-function(v){

  type=2 #type of design
  m<-((v-1)/2)
  k=0
  x1<-matrix(,nrow=1,ncol=0)
  while(k<50){
    x2<-((2*k)+1)


    k=k+1
    x1<-cbind(x1,x2)
  }
  x1<-x1[c(-1,-5,-13,-14,-25,-41)] #domain of tratments
  v= #Number of treatments

    if(any(x1==v)){
      i=1
      x<-matrix(,nrow=1, ncol=0)
      while(i<=m){
        x1<-t(c(i,-i))

        i=i+1
        x<-cbind(x,x1)
      }
      #print(x)
      x<-cbind(0,x)
      col1<-x[,-v]%%(v-1)
      mat<-matrix(0,nrow=v-1,ncol=v-1)
      mat<-rbind(col1,mat)
      k=1
      while(k<v){
        mat[(k+1),]<-(mat[k,]+1)
        xx<-mat
        k=k+1
      }

      mat<-(xx%%(v)+1)
      matt<-mat
      c1<-matrix("|",nrow=v,ncol=1)
      mat<-cbind(mat[,(v-1)],c1,mat,c1,mat[,1])
      print("NBBD with left and right border plots", quote=FALSE)
      prmatrix(mat,rowlab=rep("",v),collab=rep("",(v-1)+4), quote = FALSE)

      #break

      b<-v

      print(c("Number of blocks",b), quote = FALSE)
      #number_of_replication
      r<-v-1
      print(c("Number of replications",r), quote = FALSE)
      #number_of_treatments
      print(c("Number of treatments",v), quote = FALSE)
      # block_size
      k<-v-1
      print(c("Block size",k), quote = FALSE)
      #both_side_neighbour_appears
      u<-1
      print(c(" Number of times each treatment appearing as both left and right neighbour to each other treatments",u) ,quote= FALSE)

      print("The final designs are neighbour-balanced at distance 2 when v is prime, as shown by Lawless (1971) and Keedwell (1984).",quote = FALSE)

      ################################



      mat1<-matt
      mat2<-cbind(mat1[,ncol(mat1)],mat1[,c(1:(ncol(mat1)-1))])
      mat3<-cbind(mat1[,c(2:(ncol(mat1)),1)])
      ################################
      incident1<-matrix(0,nrow=length(mat1), ncol=v)

      i=1
      while(i<=v){
        x2<-c(which(t(mat1)==i))
        for(j in x2){
          incident1[j,i]<-(incident1[j,i]+1)
        }
        i=i+1
      }
      #
      #print(incident1)
      ################################
      incident2<-matrix(0,nrow=length(mat1), ncol=v)
      i=1
      while(i<=v){
        x2<-c(which(t(mat2)==i))
        for(j in x2){
          incident2[j,i]<-incident2[j,i]+1
        }
        i=i+1
      }
      #

      #print(incident2)
      ##############################
      incident3<-matrix(0,nrow=length(mat1), ncol=v)
      i=1
      while(i<=v){
        x2<-c(which(t(mat3)==i))
        for(j in x2){
          incident3[j,i]<-incident3[j,i]+1
        }
        i=i+1
      }
      #

      #print(incident3)
      #####################################
      #############################
      #D_matrix

      k=1
      d_mat<-matrix(,nrow=0,ncol=b)
      while(k<=b){
        xd<-matrix(,nrow=length(mat1[k,]),ncol=0)
        id<-matrix(0,nrow=length(mat1[k,]),ncol=b)
        id[,k]=1
        xd<-id
        d_mat<-rbind(d_mat,(xd))
        #print(d_mat)

        k=k+1

      }
      #d_mat
      ###################################

      ##################################
      x1_mat<-cbind(incident1,incident2,incident3)
      vec1n<-matrix(1,nrow=nrow(incident3),ncol=1)
      x2_mat<-cbind(vec1n,d_mat)
      #################x1 prime x1
      x1_prime_x1<-t(x1_mat)%*% x1_mat
      #################x1 prime x2
      x1_prime_x2<-t(x1_mat)%*% x2_mat
      #################x2 prime x2
      x2_prime_x2<-t(x2_mat)%*% x2_mat

      #########################



      #joint C matrix
      Cmatrix<-(x1_prime_x1)-(x1_prime_x2)%*%MASS::ginv(x2_prime_x2)%*%(t(x1_prime_x2))
      Cmatrix<-round(Cmatrix,digits = 3)
      #print("Joint C matrix", quote=FALSE)
      #print(Cmatrix)


      #The information matrix for estimating the direct effects of treatments
      c11<-Cmatrix[c(1:v),c(1:v)]
      c12<-Cmatrix[c(1:v),c((v+1):(3*v))]
      c21<-t(c12)
      c22<-Cmatrix[c((v+1):(3*v)),c((v+1):(3*v))]

      C_tau<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the direct effects of treatments", quote=FALSE)
      print(round(C_tau,digits = 3))
      #The information matrix for estimating the left neighbour effects of treatments

      c11<-Cmatrix[c((v+1):(2*v)),c((v+1):(2*v))]
      c121<-Cmatrix[c((v+1):(2*v)),c((1:v))]
      c122<-Cmatrix[c((v+1):(2*v)),c((2*v+1):(3*v))]
      c12<-cbind(c121,c122)
      c21<-t(c12)
      c221<-Cmatrix[c(1:v),c(1:v)]
      c222<-Cmatrix[c(1:v),c(((2*v)+1):(3*v))]
      c223<-Cmatrix[c(((2*v)+1):(3*v)),c(1:v)]
      c224<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
      c225<-cbind(c221,c222)
      c226<-cbind(c223,c224)
      c22<-rbind(c225,c226)
      C_ro<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the left neighbour effects of treatments", quote=FALSE)
      print(round(C_ro,digits = 3))
      #The information matrix for estimating the right neighbour effects of treatments
      c11<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
      c12<-Cmatrix[c(((2*v)+1):(3*v)),c(1:(2*v))]
      c21<-t(c12)
      c22<-Cmatrix[c(1:(2*v)),c(1:(2*v))]

      C_del<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the right neighbour effects of treatments", quote=FALSE)
      print(round(C_del,digits = 3))

    }else{
      print(" Plese enter a correct value",quote=FALSE)
    }
}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


#' NBB Design  of Type 3 (NBBD 3)
#'@name nbbd3
#' @param v Number of treatments (When v is even but not a power of 2)
#'@description When v is even number but not a power of 2 (v<=100), the function generates a class of neighbour balanced block designs for a given number of treatments. This class of design is generated using the method of constructing this designs which is listed as Type 2 designs of (Azais, 1993). It also gives the parameters of the design, information matrix for estimating the contrast pertaining to direct and neighbour effects (both left and right) of the treatments.
#' @return It generates neighbour balanced block designs when v is an even number but not a power of 2. For example v = 6, 10, 12, 14…
#' @export
#'@references Azais, J.M., Bailey, R.A. and Monod, H. (1993)<DOI: 10.2307/2532269>."A catalogue of efficient neighbour designs with border plots".
#' @examples
#' library(NBBDesigns)
#' nbbd3(6)
nbbd3<-function(v){

  vv=v-1
  type=2
  m<-((vv-1)/2)
  k=1
  x1<-matrix(,nrow=1,ncol=0)
  while(k<50){
    x2<-((2*k))
    #print(odd)

    k=k+1
    x1<-cbind(x1,x2)
  }


  x1<-x1[-c(2,4,8,16,32,64)]
  size<-length(x1)

  v=  #Number of treatments
    if(v %in% x1){
      flag<-TRUE

      i=1
      x<-matrix(,nrow=1, ncol=0)
      while(i<=m){
        x1<-t(c(i,-i))

        i=i+1
        x<-cbind(x,x1)
      }

      x<-cbind(0,x)
      col1<-x[,-vv]%%(vv-1)
      mat<-matrix(0,nrow=vv-1,ncol=vv-1)
      mat<-rbind(col1,mat)
      k=1
      while(k<vv){
        mat[(k+1),]<-(mat[k,]+1)
        xx<-mat
        k=k+1
      }
      mat<-(xx%%(vv))+1
      matt<-mat
      c2<-matrix(v,nrow=vv,ncol=1)
      c1<-matrix("|",nrow=vv,ncol=1)
      gl<-c(1:vv)
      mat<-cbind(mat[,(vv-1)],c1,mat,c1,mat[,1])
      mat1<-mat[,c(1,2,3)]
      mat2<-mat[,-c(1,2,3)]

      mat1<-cbind(mat1,c2)
      mat<-cbind(mat1,mat2)
      mat3<-mat[,-c(1,2,(v+2),(v+3))]
      mat3<-rbind(mat3,gl)
      c1<-matrix("|",nrow=v,ncol=1)
      mat3<-cbind(mat3[,vv],c1,mat3,c1,mat3[,1])

      print("NBBD with left and right border plots", quote=FALSE)
      prmatrix(mat3,rowlab=rep("",v),collab=rep("",(v-1)+4), quote = FALSE)


      b<-v


      print(c("Number of blocks",b), quote = FALSE)
      #number_of_replication
      r<-v-1
      print(c("Number of replications",r), quote = FALSE)
      #number_of_treatments
      print(c("Number of treatments",v), quote = FALSE)
      # block_size
      k<-v-1
      print(c("Block size",k), quote = FALSE)
      #both_side_neighbour_appears
      u<-1
      print(c(" Number of times each treatment appearing as both left and right neighbour to each other treatments",u) ,quote= FALSE)




      mat1<-matt
      mat2<-cbind(mat1[,ncol(mat1)],mat1[,c(1:(ncol(mat1)-1))])
      mat3<-cbind(mat1[,c(2:(ncol(mat1)),1)])
      ################################
      incident1<-matrix(0,nrow=length(mat1), ncol=vv)

      i=1
      while(i<=vv){
        x2<-c(which(t(mat1)==i))
        for(j in x2){
          incident1[j,i]<-(incident1[j,i]+1)
        }
        i=i+1
      }
      #

      ################################
      incident2<-matrix(0,nrow=length(mat1), ncol=vv)
      i=1
      while(i<=vv){
        x2<-c(which(t(mat2)==i))
        for(j in x2){
          incident2[j,i]<-incident2[j,i]+1
        }
        i=i+1
      }
      #

      ##############################
      incident3<-matrix(0,nrow=length(mat1), ncol=vv)
      i=1
      while(i<=vv){
        x2<-c(which(t(mat3)==i))
        for(j in x2){
          incident3[j,i]<-incident3[j,i]+1
        }
        i=i+1
      }
      #

      #####################################
      #############################
      #D_matrix

      k=1
      d_mat<-matrix(,nrow=0,ncol=b)
      while(k<=(v-1)){
        xd<-matrix(,nrow=length(mat1[k,]),ncol=0)
        id<-matrix(0,nrow=length(mat1[k,]),ncol=b)
        id[,k]=1
        xd<-id
        d_mat<-rbind(d_mat,(xd))
        #print(d_mat)

        k=k+1

      }

      ###################################

      ##################################
      x1_mat<-cbind(incident1,incident2,incident3)
      vec1n<-matrix(1,nrow=nrow(incident3),ncol=1)
      x2_mat<-cbind(vec1n,d_mat)
      #################x1 prime x1
      x1_prime_x1<-t(x1_mat)%*% x1_mat
      #################x1 prime x2
      x1_prime_x2<-t(x1_mat)%*% x2_mat
      #################x2 prime x2
      x2_prime_x2<-t(x2_mat)%*% x2_mat

      #########################

      #########################G inverse


      #joint C matrix

      Cmatrix<-(x1_prime_x1)-(x1_prime_x2)%*%MASS::ginv(x2_prime_x2)%*%(t(x1_prime_x2))
      Cmatrix<-round(Cmatrix,digits = 3)
     # print("Joint C matrix", quote=FALSE)
      #print(Cmatrix)

      ########################

      #The information matrix for estimating the direct effects of treatments
      c11<-Cmatrix[c(1:vv),c(1:vv)]
      c12<-Cmatrix[c(1:vv),c((vv+1):(3*vv))]
      c21<-t(c12)
      c22<-Cmatrix[c((vv+1):(3*vv)),c((vv+1):(3*vv))]

      C_tau<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the direct effects of treatments", quote=FALSE)
      #print(C_tau)
      print(round(C_tau,digits = 3))
      #The information matrix for estimating the left neighbour effects of treatments

      c11<-Cmatrix[c((vv+1):(2*vv)),c((vv+1):(2*vv))]
      c121<-Cmatrix[c((vv+1):(2*vv)),c((1:vv))]
      c122<-Cmatrix[c((vv+1):(2*vv)),c((2*vv+1):(3*vv))]
      c12<-cbind(c121,c122)
      c21<-t(c12)
      c221<-Cmatrix[c(1:vv),c(1:vv)]
      c222<-Cmatrix[c(1:vv),c(((2*vv)+1):(3*vv))]
      c223<-Cmatrix[c(((2*vv)+1):(3*vv)),c(1:vv)]
      c224<-Cmatrix[c(((2*vv)+1):(3*vv)),c(((2*vv)+1):(3*vv))]
      c225<-cbind(c221,c222)
      c226<-cbind(c223,c224)
      c22<-rbind(c225,c226)
      C_ro<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the contrast pertaining to the left neighbour effects of treatments", quote=FALSE)
      #print(C_ro)
      print(round(C_ro,digits = 3))
      #The information matrix for estimating the right neighbour effects of treatments
      c11<-Cmatrix[c(((2*vv)+1):(3*vv)),c(((2*vv)+1):(3*vv))]
      c12<-Cmatrix[c(((2*vv)+1):(3*vv)),c(1:(2*vv))]
      c21<-t(c12)
      c22<-Cmatrix[c(1:(2*vv)),c(1:(2*vv))]
      print("The information matrix for estimating the contrast pertaining to the contrast pertaining to the right neighbour effects of treatments", quote=FALSE)
      C_del<-c11-c12%*%MASS::ginv(c22)%*%c21
     # print(C_del)
      print(round(C_del,digits = 3))
      #######################

    }else{
      print("Please enter a correct value")
    }
}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


#' PNBB Design of Type 1 (PNBBD 1)
#'@name pnbbd1
#' @param v Number of treatments (v), v>4 should be of the form v=s^2. For example v= 8, 9…..
#'@description
#'A block design with neighbour effects is said to be partially neighbour balanced based on m-class association scheme if two treatments 'Theta' and 'Phi' that are mutually u-th associates (u = 1, 2,…, m) appear as neighbours (left and right) 'Mu'_1u times. For this class of design 'Mu'_11=1 and 'Mu'_12=0. The number of first associates is (v/2) and number of second associates is (v/2 -1).
#'When v is a power of 2, the function will generate a class of partially neighbour balanced block designs. It also gives the parameters of the design, information matrix for estimating the contrast pertaining to direct and neighbour effects (both left and right) of the treatments.
#' @export
#' @references Azais, J.M., Bailey, R.A. and Monod, H. (1993)<DOI: 10.2307/2532269>."A catalogue of efficient neighbour designs with border plots".
#'@note Here v should be greater than 4 i.e v>4.
#' @examples
#' library(NBBDesigns)
#' pnbbd1(8)
pnbbd1<-function(v){

  p=((log(v))/log(2))
  type=1
  v=
    if(v%%4==0 && v>4){
      k=((2^(p-1))-1)
      x1<-matrix(,nrow=1,ncol=0)
      while(k>=0){
        x2<-((2*k)+1)
        #print(odd)

        k=k-1
        x1<-cbind(x1,x2)
      }
      #x1

      mat<-matrix(,nrow=0,ncol=v)
      for(i in x1){
        x2<-matrix(,nrow=1,ncol=0)
        x3<-t(seq(0,((i*v)-1), by=i))
        x2<-cbind(x2,x3)
        mat<-(rbind((x2),(mat)))
        i=i+1
      }
      mat<-(mat%%v+1)
      matt<-mat
      c1<-matrix("|",nrow=length(x1),ncol=1)
      mat<-cbind(mat[,(v)],c1,mat,c1,mat[,1])

      print("PNBBD with left and right border plots", quote=FALSE)
      prmatrix(mat,rowlab=rep("",v),collab=rep("",(v)+4), quote = FALSE)

      b<-(2^(p-1))

      print(c("Number of blocks",b), quote = FALSE)
      #number_of_replication
      r<-(2^(p-1))
      print(c("Number of replications",r), quote = FALSE)
      #number_of_treatments
      print(c("Number of treatments",v), quote = FALSE)
      # block_size
      k<-v
      print(c("Block size",k), quote = FALSE)
      #both_side_neighbour_appears
      u<-1
      print(c(" Number of times each treatment appearing as both left and right neighbour to each other treatments",u) ,quote= FALSE)
      print(c(" Number of first associates",(v/2)) ,quote= FALSE)
      print(c(" Number of second associates",(v/2)-1) ,quote= FALSE)


      ########################

      mat1<-matt
      mat2<-cbind(mat1[,ncol(mat1)],mat1[,c(1:(ncol(mat1)-1))])
      mat3<-cbind(mat1[,c(2:(ncol(mat1)),1)])
      ################################
      incident1<-matrix(0,nrow=length(mat1), ncol=v)

      i=1
      while(i<=v){
        x2<-c(which(t(mat1)==i))
        for(j in x2){
          incident1[j,i]<-(incident1[j,i]+1)
        }
        i=i+1
      }
      #
      #print(incident1)
      ################################
      incident2<-matrix(0,nrow=length(mat1), ncol=v)
      i=1
      while(i<=v){
        x2<-c(which(t(mat2)==i))
        for(j in x2){
          incident2[j,i]<-incident2[j,i]+1
        }
        i=i+1
      }
      #

      #print(incident2)
      ##############################
      incident3<-matrix(0,nrow=length(mat1), ncol=v)
      i=1
      while(i<=v){
        x2<-c(which(t(mat3)==i))
        for(j in x2){
          incident3[j,i]<-incident3[j,i]+1
        }
        i=i+1
      }
      #

      #print(incident3)
      #####################################
      #############################
      #D_matrix

      k=1
      d_mat<-matrix(,nrow=0,ncol=b)
      while(k<=b){
        xd<-matrix(,nrow=length(mat1[k,]),ncol=0)
        id<-matrix(0,nrow=length(mat1[k,]),ncol=b)
        id[,k]=1
        xd<-id
        d_mat<-rbind(d_mat,(xd))
        #print(d_mat)

        k=k+1

      }
      #d_mat
      ###################################

      ##################################
      x1_mat<-cbind(incident1,incident2,incident3)
      vec1n<-matrix(1,nrow=nrow(incident3),ncol=1)
      x2_mat<-cbind(vec1n,d_mat)
      #################x1 prime x1
      x1_prime_x1<-t(x1_mat)%*% x1_mat
      #################x1 prime x2
      x1_prime_x2<-t(x1_mat)%*% x2_mat
      #################x2 prime x2
      x2_prime_x2<-t(x2_mat)%*% x2_mat

      #########################



      #x2_prime_x2%*%ginv%*%x2_prime_x2
      #joint C matrix
      Cmatrix<-(x1_prime_x1)-(x1_prime_x2)%*%MASS::ginv(x2_prime_x2)%*%(t(x1_prime_x2))
      Cmatrix<-round(Cmatrix,digits = 3)
      #print("Joint C matrix", quote=FALSE)
      #print(Cmatrix)


      #The information matrix for estimating the direct effects of treatments
      c11<-Cmatrix[c(1:v),c(1:v)]
      c12<-Cmatrix[c(1:v),c((v+1):(3*v))]
      c21<-t(c12)
      c22<-Cmatrix[c((v+1):(3*v)),c((v+1):(3*v))]

      C_tau<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the direct effects of treatments", quote=FALSE)
      print(round(C_tau,digits = 3))
      #The information matrix for estimating the left neighbour effects of treatments

      c11<-Cmatrix[c((v+1):(2*v)),c((v+1):(2*v))]
      c121<-Cmatrix[c((v+1):(2*v)),c((1:v))]
      c122<-Cmatrix[c((v+1):(2*v)),c((2*v+1):(3*v))]
      c12<-cbind(c121,c122)
      c21<-t(c12)
      c221<-Cmatrix[c(1:v),c(1:v)]
      c222<-Cmatrix[c(1:v),c(((2*v)+1):(3*v))]
      c223<-Cmatrix[c(((2*v)+1):(3*v)),c(1:v)]
      c224<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
      c225<-cbind(c221,c222)
      c226<-cbind(c223,c224)
      c22<-rbind(c225,c226)
      C_ro<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the left neighbour effects of treatments", quote=FALSE)
      print(round(C_ro,digits = 3))
      #The information matrix for estimating the right neighbour effects of treatments
      c11<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
      c12<-Cmatrix[c(((2*v)+1):(3*v)),c(1:(2*v))]
      c21<-t(c12)
      c22<-Cmatrix[c(1:(2*v)),c(1:(2*v))]

      C_del<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the right neighbour effects of treatments", quote=FALSE)
      print(round(C_del,digits = 3))


    }else{
      print(" Plese enter a correct value",quote=FALSE)
    }
}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




#' PNBB Design of Type 2 (PNBBD 2)
#'@name pnbbd2
#' @param v Number of treatments (v), v should be a prime number
#'@description
#'A block design with neighbour effects is said to be partially neighbour balanced based on m-class association scheme if two treatments 'Theta' and 'Phi' that are mutually u-th associates (u = 1, 2,…, m) appear as neighbours (left and right) 'Mu'_1u times. The design so obtained is a ((v-1)/2) associate classes partially variance balanced design following a varying circular association scheme. When v is a prime number, the function will generate a class of partially neighbour balanced block designs. It also gives the parameters of the design, information matrix for estimating the contrast pertaining to direct and neighbour effects (both left and right) of the treatments.
#' @return It gives Partially Neighbour Balanced Block Designs for v, when v is any prime number.
#' @export
#' @references Azais, J.M., Bailey, R.A. and Monod, H. (1993)<DOI: 10.2307/2532269>."A catalogue of efficient neighbour designs with border plots".
#'@note v should be greater than 5 i.e v>5.
#' @examples
#' library(NBBDesigns)
#' pnbbd2(7)
pnbbd2<-function(v){
  p=(v-1)/2
  type=2

  i=2
  while(i<=(v/2)){
    if(v%%i!=0){
      i=i+1
    }else{
      #print(c(v, "Entered number is not a prime "), quote=FALSE)
      break
    }
  }
  v=
    if(i>v/2 && v>5){
      #print(c(v, "is prime number"), quote=FALSE)

      x1<-c(p:1)

      mat<-matrix(,nrow=0,ncol=v)
      for(i in x1){
        x2<-matrix(,nrow=1,ncol=0)
        x3<-t(seq(0,((i*v)-1), by=i))
        x2<-cbind(x2,x3)
        mat<-(rbind((x2),(mat)))
        i=i+1
      }
      mat<-(mat%%v+1)
      matt<-mat
      c1<-matrix("|",nrow=length(x1),ncol=1)
      mat<-cbind(mat[,(v)],c1,mat,c1,mat[,1])

      print("PNBBD with left and right border plots", quote=FALSE)
      prmatrix(mat,rowlab=rep("",v),collab=rep("",(v)+4), quote = FALSE)
      b<-p

      print(c("Number of blocks",b), quote = FALSE)
      #number_of_replication
      r<-p
      print(c("Number of replications",r), quote = FALSE)
      #number_of_treatments
      print(c("Number of treatments",v), quote = FALSE)
      # block_size
      k<-v
      print(c("Block size",k), quote = FALSE)
      #both_side_neighbour_appears
      # u<-1
      # print(c(" Number of times each treatment appearing as both left and right neighbour to each other treatments",u) ,quote= FALSE)
      #

      ###########################

      mat1<-matt
      mat2<-cbind(mat1[,ncol(mat1)],mat1[,c(1:(ncol(mat1)-1))])
      mat3<-cbind(mat1[,c(2:(ncol(mat1)),1)])
      ################################
      incident1<-matrix(0,nrow=length(mat1), ncol=v)

      i=1
      while(i<=v){
        x2<-c(which(t(mat1)==i))
        for(j in x2){
          incident1[j,i]<-(incident1[j,i]+1)
        }
        i=i+1
      }
      #
      #print(incident1)
      ################################
      incident2<-matrix(0,nrow=length(mat1), ncol=v)
      i=1
      while(i<=v){
        x2<-c(which(t(mat2)==i))
        for(j in x2){
          incident2[j,i]<-incident2[j,i]+1
        }
        i=i+1
      }
      #

      #print(incident2)
      ##############################
      incident3<-matrix(0,nrow=length(mat1), ncol=v)
      i=1
      while(i<=v){
        x2<-c(which(t(mat3)==i))
        for(j in x2){
          incident3[j,i]<-incident3[j,i]+1
        }
        i=i+1
      }
      #

      #print(incident3)
      #####################################
      #############################
      #D_matrix

      k=1
      d_mat<-matrix(,nrow=0,ncol=b)
      while(k<=b){
        xd<-matrix(,nrow=length(mat1[k,]),ncol=0)
        id<-matrix(0,nrow=length(mat1[k,]),ncol=b)
        id[,k]=1
        xd<-id
        d_mat<-rbind(d_mat,(xd))
        #print(d_mat)

        k=k+1

      }
      #d_mat
      ###################################

      ##################################
      x1_mat<-cbind(incident1,incident2,incident3)
      vec1n<-matrix(1,nrow=nrow(incident3),ncol=1)
      x2_mat<-cbind(vec1n,d_mat)
      #################x1 prime x1
      x1_prime_x1<-t(x1_mat)%*% x1_mat
      #################x1 prime x2
      x1_prime_x2<-t(x1_mat)%*% x2_mat
      #################x2 prime x2
      x2_prime_x2<-t(x2_mat)%*% x2_mat

      #########################



      #x2_prime_x2%*%ginv%*%x2_prime_x2
      #joint C matrix
      Cmatrix<-(x1_prime_x1)-(x1_prime_x2)%*%MASS::ginv(x2_prime_x2)%*%(t(x1_prime_x2))
      Cmatrix<-round(Cmatrix,digits = 3)
     # print("Joint C matrix", quote=FALSE)
      #print(Cmatrix)


      #The information matrix for estimating the direct effects of treatments
      c11<-Cmatrix[c(1:v),c(1:v)]
      c12<-Cmatrix[c(1:v),c((v+1):(3*v))]
      c21<-t(c12)
      c22<-Cmatrix[c((v+1):(3*v)),c((v+1):(3*v))]

      C_tau<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the direct effects of treatments", quote=FALSE)
      print(round(C_tau,digits = 3))
      #The information matrix for estimating the left neighbour effects of treatments

      c11<-Cmatrix[c((v+1):(2*v)),c((v+1):(2*v))]
      c121<-Cmatrix[c((v+1):(2*v)),c((1:v))]
      c122<-Cmatrix[c((v+1):(2*v)),c((2*v+1):(3*v))]
      c12<-cbind(c121,c122)
      c21<-t(c12)
      c221<-Cmatrix[c(1:v),c(1:v)]
      c222<-Cmatrix[c(1:v),c(((2*v)+1):(3*v))]
      c223<-Cmatrix[c(((2*v)+1):(3*v)),c(1:v)]
      c224<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
      c225<-cbind(c221,c222)
      c226<-cbind(c223,c224)
      c22<-rbind(c225,c226)
      C_ro<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the left neighbour effects of treatments", quote=FALSE)
      print(round(C_ro,digits = 3))
      #The information matrix for estimating the right neighbour effects of treatments
      c11<-Cmatrix[c(((2*v)+1):(3*v)),c(((2*v)+1):(3*v))]
      c12<-Cmatrix[c(((2*v)+1):(3*v)),c(1:(2*v))]
      c21<-t(c12)
      c22<-Cmatrix[c(1:(2*v)),c(1:(2*v))]

      C_del<-c11-c12%*%MASS::ginv(c22)%*%c21
      print("The information matrix for estimating the contrast pertaining to the right neighbour effects of treatments", quote=FALSE)
      print(round(C_del,digits = 3))

    }else{
      print("Please enter a correct value",quote=FALSE)
    }
}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
globalVariables(c("drop1","lm","a"))
##Analysis
#' Analysis of  data
#'
#' @param data The data file should be in csv format. The columns should be named as block,treatment, left_neighbour, right_neighbour and yield as given in the example data set.
#'@description This function provides the Analysis of Variance (Type III) of the data generated from experiments conducted using a neighbour balanced/partially neighbour balanced block design.
#' @return It provides the ANOVA table.
#' @export
#'@importFrom utils read.csv
#'@importFrom stats anova
#'@examples
#'\dontrun{
#'library(NBBDesigns)
#'data<-file.choose()
#'data<-read.csv(data,header=TRUE,colClasses = c("factor","factor","factor","factor","numeric"))
#'fix(data)
#'anlys(data)
#'}
anlys<-function(data){
  model1<-lm(yield~block+treatment+left_neighbour+right_neighbour, data=data)
  a<-drop1(model1,~.,test = "F")
  aa<-a[-1,-c(3,4)]
  c<-anova(model1)
  b<-c[5,]
  b<-b[,-3]
  b<-data.frame(b)
  colnames(b)<-NULL
  colnames(b)<-colnames(aa)
  rbind.data.frame(aa,b)
}


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



