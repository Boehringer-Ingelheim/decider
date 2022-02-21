#'@keywords internal
plot_decisions_jointBLRM_int <- function(
  summary,
  type,
  probs,
  esc.rule,
  dosing.intervals,
  loss.weights,
  dloss.weights,
  ref.probs.all,
  ewoc.threshold,
  combi.heatmap,
  p.int.probs,
  file.name,
  path,
  file.format,
  width=NULL,
  height=NULL,
  unit=NULL,
  underint.deact = FALSE
){
  summary_int <- summary[, (2+length(probs)+1):(2+length(probs)+length(dosing.intervals)+1)]
  if(tolower(esc.rule)%in%c("ewoc")){
    #perform ewoc escalation
    plots <- plot_decisions_jointBLRM_ewoc(summary=summary_int,
                                          type=type,
                                          heatmap=combi.heatmap,
                                          ewoc.threshold=ewoc.threshold,
                                          file.name=file.name,
                                          path=path,
                                          file.format=file.format,
                                          underint.deact = underint.deact,
                                          width=width,
                                          height=height,
                                          unit=unit)
    return(plots)

  }else if(tolower(esc.rule)%in%c("loss")){
     #perform usual loss escalation
     plots <- plot_decisions_jointBLRM_loss(
                   summary=summary_int,
                   type=type,
                   heatmap=combi.heatmap,
                   loss.weights=loss.weights,
                   p.int.probs = p.int.probs,
                   file.name=file.name,
                   path=path,
                   file.format=file.format,
                   width=width,
                   height=height,
                   unit=unit)
    return(plots)
    #return(NULL)

  }else if(tolower(esc.rule)%in%c("dynamic.loss")){
    #perform dynamic loss escalation
    plots <- plot_decisions_jointBLRM_dloss(
      summary=summary_int,
      type=type,
      probs.ref = ref.probs.all,
      heatmap=combi.heatmap,
      dloss.weights=dloss.weights,
      p.int.probs = p.int.probs,
      file.name=file.name,
      path=path,
      file.format=file.format,
      width=width,
      height=height,
      unit=unit)
    return(plots)
  }
}



#----------------------------------------------------
#----------------------------------------------------
#EWOC escalation
#----------------------------------------------------
#----------------------------------------------------


#----------------------------------------------------
#general plotting function for EWOC escalation
#----------------------------------------------------

#'@keywords internal
plot_decisions_jointBLRM_ewoc <- function(
  summary,
  type,
  heatmap,
  ewoc.threshold,
  file.name,
  path,
  file.format,
  underint.deact = FALSE,
  width=NULL,
  height=NULL,
  unit=NULL
){
  if(type %in% c("mono1", "mono2", "combi")){
    return(plot_decisions_jointBLRM_ewoct(summary=summary, type=type, heatmap=heatmap, ewoc.threshold=ewoc.threshold,
                                          file.name=file.name,
                                          path=path,
                                          file.format=file.format,
                                          underint.deact = underint.deact,
                                          width=width,
                                          height=height,
                                          unit=unit))
  }else{
    #plot type is "all"

    #Determine if mono1, mono2 and combi doses are in summary,
    #and call the corresponding type-specific plot functions for each
    #detected type.

    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose1 <- rep(0, times=ndoses)
    dose2 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)

    #count number of dose levels per type
    nmon1 <- 0
    nmon2 <- 0
    ncom <- 0
    mono1 <- FALSE
    mono2 <- FALSE
    combi <- FALSE
    for(d in 1:ndoses){
      dose1[d] <- as.numeric(strdose[[d]][1])
      dose2[d] <- as.numeric(strdose[[d]][2])
      if(dose1[d]>0 & dose2[d]==0){
        nmon1 <- nmon1+1
        mono1<-TRUE
      }else if(dose1[d]==0 & dose2[d]>0){
        nmon2 <- nmon2+1
        mono2<-TRUE
      }else if(dose1[d]>0 & dose2[d]>0){
        ncom <- nmon2+1
        combi<-TRUE
      }
    }

    nacc <- mono1+mono2+combi
    plots <- list()

    if(mono1){
      summ <- summary[which(dose1>0 & dose2==0), ]
      plmon1 <- plot_decisions_jointBLRM_ewoct(summary=summ,
                                               type="mono1",
                                               heatmap = heatmap,
                                               ewoc.threshold=ewoc.threshold,
                                               file.name=file.name,
                                               path=path,
                                               file.format=file.format,
                                               underint.deact = underint.deact,
                                               width=width,
                                               height=height,
                                               unit=unit)
      if(nacc==1){
        return(plmon1)
      }else{
        plots[["Mono1"]] <- plmon1
      }
    }

    if(mono2){
      summ <- summary[which(dose2>0 & dose1==0), ]
      plmon2 <- plot_decisions_jointBLRM_ewoct(summary=summ,
                                               type="mono2",
                                               heatmap = heatmap,
                                               ewoc.threshold=ewoc.threshold,
                                               file.name=file.name,
                                               path=path,
                                               file.format=file.format,
                                               underint.deact = underint.deact,
                                               width=width,
                                               height=height,
                                               unit=unit)
      if(nacc==1){
        return(plmon2)
      }else{
        plots[["Mono2"]] <- plmon2
      }
    }

    if(combi){
      summ <- summary[which(dose2>0 & dose1>0), ]
      plcom <- plot_decisions_jointBLRM_ewoct(summary=summ,
                                               type="combi",
                                               heatmap = heatmap,
                                              ewoc.threshold=ewoc.threshold,
                                              file.name=file.name,
                                              path=path,
                                              file.format=file.format,
                                              underint.deact = underint.deact,
                                              width=width,
                                              height=height,
                                              unit=unit)
      if(nacc==1){
        return(plcom)
      }else{
        plots[["Combi"]] <- plcom
      }
    }

    return(plots)

  }

}

#----------------------------------------------------
#Plotting function for EWOC escalation with types
#in "mono1", "mono2" or "combi".
#different function used for "all".
#----------------------------------------------------

#'@keywords internal
plot_decisions_jointBLRM_ewoct <- function(
  summary,
  type,
  heatmap,
  ewoc.threshold,
  file.name,
  path,
  file.format,
  underint.deact = FALSE,
  width=NULL,
  height=NULL,
  unit=NULL
){

  if(file.format=="pdf"){
    ffint <- cairo_pdf
  }else{
    ffint <- file.format
  }
  if(type=="mono1"){

    #----------------------------------------------------
    #Plots for EWOC escalation: Mono 1
    #----------------------------------------------------
    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose1 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)
    for(d in 1:ndoses){
      dose1[d] <- as.numeric(strdose[[d]][1])
    }

    #create data frame for plotting
    ldat <- list()
    ewoc <- rep(0, times = ndoses)
    ewoc.sat <- FALSE
    ewoc.not.sat <- FALSE
    for(i in 1:ndoses){
      if(summary[i,3]>=ewoc.threshold){
        ewoc[i]<- paste0("P(Over)\u2265", ewoc.threshold)
        ewoc.not.sat <- TRUE
      }else{
        ewoc[i]<- paste0("P(Over)<", ewoc.threshold)
        ewoc.sat <- TRUE
      }
    }
    ldat[["dose"]] <- factor(dose1, levels=unique(dose1))
    ldat[["under"]] <- summary[,1]
    ldat[["target"]] <- summary[, 2]
    ldat[["over"]] <- summary[, 3]
    ldat[["ewoc"]]<- ewoc
    datf <- as.data.frame(ldat)

    if(ewoc.not.sat){
      if(ewoc.sat){
        color_vec <- c(  "tomato1","firebrick", "forestgreen","gold")
      }else{
        color_vec <- c( "firebrick", "forestgreen","gold")
      }
    }else{
      color_vec <- c(  "tomato1", "forestgreen","gold")
    }

    if(!underint.deact){
      punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
        geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                            size = 2, linetype = "solid"),
            panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                            colour = "white"),
            panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                            colour = "white"),
            plot.margin = unit(c(-1, 3, 3, 3), "mm"),
            plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
            plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
            axis.title = element_text(size=12, face = "bold", family= "sans"),
            axis.text = element_text(size=10, family= "sans"),
            legend.key.size = unit(5, "mm"),
            legend.title = element_text(size=12, family= "sans"),
            legend.text = element_text(size=10, family= "sans"),
            legend.position = "top")
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill = "ewoc"))+
        geom_bar(stat="identity", aes(fill = "P(Under)"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec)+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        geom_hline(aes(yintercept=ewoc.threshold, linetype ="EWOC"), color = "firebrick", size = 0.8)+
        scale_linetype_manual("", values = "dotted")+
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
              #            colour = "white"),
              plot.margin = unit(c(3, 3, -2, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.margin = margin(0, 0, -3, 0, "mm"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
              legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                              colour = "white"),
              plot.margin = unit(c(-1, 3, -1, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans"),
              legend.position = "top")

      plist <- list(pover, ptar, punder)
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_mono1.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=3),
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_mono1.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=3),
              dpi="retina",
              width=width[1],
              height=height[1],
              units = unit
            )

          }
        }
      }
      return(grid.arrange(grobs=plist, nrow=3))

    }else{
      #Underdosing is deactivated.
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill = "ewoc"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec[1:(length(color_vec)-1)])+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        geom_hline(aes(yintercept=ewoc.threshold, linetype ="EWOC"), color = "firebrick", size = 0.8)+
        scale_linetype_manual("", values = "dotted")+
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
              #            colour = "white"),
              plot.margin = unit(c(3, 3, -2, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.margin = margin(0, 0, -3, 0, "mm"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
              legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                              colour = "white"),
              plot.margin = unit(c(-1, 3, 3, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans"),
              legend.position = "top")

      plist <- list(pover, ptar)
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_mono1.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=2),
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_mono1.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=2),
              dpi="retina",
              width=width[1],
              height=height[1],
              units = unit
            )

          }
        }
      }
      return(grid.arrange(grobs=plist, nrow=2))


    }
  }else if(type=="mono2"){

    #----------------------------------------------------
    #Plots for EWOC escalation: Mono 2
    #----------------------------------------------------
    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose2 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)
    for(d in 1:ndoses){
      dose2[d] <- as.numeric(strdose[[d]][2])
    }

    #create data frame for plotting
    ldat <- list()
    ewoc <- rep(0, times = ndoses)
    ewoc.sat <- FALSE
    ewoc.not.sat <- FALSE
    for(i in 1:ndoses){
      if(summary[i,3]>=ewoc.threshold){
        ewoc[i]<- paste0("P(Over)\u2265", ewoc.threshold)
        ewoc.not.sat <- TRUE
      }else{
        ewoc[i]<- paste0("P(Over)<", ewoc.threshold)
        ewoc.sat <- TRUE
      }
    }
    ldat[["dose"]] <- factor(dose2, levels=unique(dose2))
    ldat[["under"]] <- summary[,1]
    ldat[["target"]] <- summary[, 2]
    ldat[["over"]] <- summary[, 3]
    ldat[["ewoc"]]<- ewoc
    datf <- as.data.frame(ldat)

    if(ewoc.not.sat){
      if(ewoc.sat){
        color_vec <- c(  "tomato1","firebrick", "forestgreen","gold")
      }else{
        color_vec <- c( "firebrick", "forestgreen","gold")
      }
    }else{
      color_vec <- c(  "tomato1", "forestgreen","gold")
    }

    if(!underint.deact){
      punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
        geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                              colour = "white"),
              plot.margin = unit(c(-1, 3, 3, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans"),
              legend.position = "top")
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill = "ewoc"))+
        geom_bar(stat="identity", aes(fill = "P(Under)"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec)+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        geom_hline(aes(yintercept=ewoc.threshold, linetype ="EWOC"), color = "firebrick", size = 0.8)+
        scale_linetype_manual("", values = "dotted")+
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
              #            colour = "white"),
              plot.margin = unit(c(3, 3, -2, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.margin = margin(0, 0, -3, 0, "mm"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
              legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                              colour = "white"),
              plot.margin = unit(c(-1, 3, -1, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans"),
              legend.position = "top")

      plist <- list(pover, ptar, punder)

      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_mono2.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=3),
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_mono2.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=3),
              dpi="retina",
              width=width[2],
              height=height[2],
              units = unit
            )
          }
        }
      }
      return(grid.arrange(grobs=plist, nrow=3))
    }else{
      #underdosing interval deactivated
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill = "ewoc"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec[1:(length(color_vec)-1)])+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        geom_hline(aes(yintercept=ewoc.threshold, linetype ="EWOC"), color = "firebrick", size = 0.8)+
        scale_linetype_manual("", values = "dotted")+
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
              #            colour = "white"),
              plot.margin = unit(c(3, 3, -2, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.margin = margin(0, 0, -3, 0, "mm"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
              legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                              size = 2, linetype = "solid"),
              panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                              colour = "white"),
              panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                              colour = "white"),
              plot.margin = unit(c(-1, 3, 3, 3), "mm"),
              plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
              plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
              axis.title = element_text(size=12, face = "bold", family= "sans"),
              axis.text = element_text(size=10, family= "sans"),
              legend.key.size = unit(5, "mm"),
              legend.title = element_text(size=12, family= "sans"),
              legend.text = element_text(size=10, family= "sans"),
              legend.position = "top")

      plist <- list(pover, ptar)

      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_mono2.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=2),
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_mono2.", file.format),
              path=path,
              device=ffint,
              plot=grid.arrange(grobs=plist, nrow=2),
              dpi="retina",
              width=width[2],
              height=height[2],
              units = unit
            )
          }
        }
      }
      return(grid.arrange(grobs=plist, nrow=2))


    }

  }else if(type=="combi"){

    #----------------------------------------------------
    #Plots for EWOC escalation: Combination
    #----------------------------------------------------

    if(!heatmap){
      #----------------------------------------------------
      #Combination plot without heatmap
      #----------------------------------------------------
      doses <- rownames(summary)
      ndoses <- length(doses)

      #create data frame for plotting
      ldat <- list()
      ewoc <- rep(0, times = ndoses)
      ewoc.sat <- FALSE
      ewoc.not.sat <- FALSE
      for(i in 1:ndoses){
        if(summary[i,3]>=ewoc.threshold){
          ewoc[i]<- paste0("P(Over)\u2265", ewoc.threshold)
          ewoc.not.sat <- TRUE
        }else{
          ewoc[i]<- paste0("P(Over)<", ewoc.threshold)
          ewoc.sat <- TRUE
        }
      }
      ldat[["dose"]] <- factor(doses, levels=unique(doses))
      ldat[["under"]] <- summary[,1]
      ldat[["target"]] <- summary[, 2]
      ldat[["over"]] <- summary[, 3]
      ldat[["ewoc"]]<- ewoc
      datf <- as.data.frame(ldat)

      if(ewoc.not.sat){
        if(ewoc.sat){
          color_vec <- c(  "tomato1","firebrick", "forestgreen","gold")
        }else{
          color_vec <- c( "firebrick", "forestgreen","gold")
        }
      }else{
        color_vec <- c(  "tomato1", "forestgreen","gold")
      }
      if(!underint.deact){
        punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
          geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="Dose Combination", y = paste0(""), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                size = 2, linetype = "solid"),
                panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                colour = "white"),
                panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                colour = "white"),
                plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                axis.title = element_text(size=12, face = "bold", family= "sans"),
                axis.text = element_text(size=10, family= "sans"),
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size=12, family= "sans"),
                legend.text = element_text(size=10, family= "sans"),
                legend.position = "top")
        pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill = "ewoc"))+
          geom_bar(stat="identity", aes(fill = "P(Under)"))+
          geom_bar(stat="identity", aes(fill = "P(Target)"))+
          geom_bar(stat="identity")+
          #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          scale_fill_manual("", values = color_vec)+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="", y = paste0(""), cex = 1.5) +
          geom_hline(aes(yintercept=ewoc.threshold, linetype ="EWOC"), color = "firebrick", size = 0.8)+
          scale_linetype_manual("", values = "dotted")+
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                size = 2, linetype = "solid"),
                panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                colour = "white"),
                panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                #            colour = "white"),
                plot.margin = unit(c(3, 3, -2, 3), "mm"),
                plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                axis.title = element_text(size=12, face = "bold", family= "sans"),
                axis.text = element_text(size=10, family= "sans"),
                legend.margin = margin(0, 0, -3, 0, "mm"),
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size=12, family= "sans"),
                legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                legend.position = "top")
        #pover
        ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
          geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="", y = paste0("Probability"), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                size = 2, linetype = "solid"),
                panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                colour = "white"),
                panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                colour = "white"),
                plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                axis.title = element_text(size=12, face = "bold", family= "sans"),
                axis.text = element_text(size=10, family= "sans"),
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size=12, family= "sans"),
                legend.text = element_text(size=10, family= "sans"),
                legend.position = "top")

        plist <- list(pover, ptar, punder)


        if(!is.null(file.name) & !is.null(path)){
          if(dir.exists(file.path(path))){
            if(is.null(unit)){
              ggsave(
                filename=paste0(file.name, "_combi.", file.format),
                path=path,
                device=ffint,
                plot=grid.arrange(grobs=plist, nrow=3),
                dpi="retina"
              )
            }else{
              ggsave(
                filename=paste0(file.name, "_combi.", file.format),
                path=path,
                device=ffint,
                plot=grid.arrange(grobs=plist, nrow=3),
                dpi="retina",
                width=width[3],
                height=height[3],
                units = unit
              )
            }
          }
        }
        return(grid.arrange(grobs=plist, nrow=3))
      }else{
        #Underdosing deactivaTED
        pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill = "ewoc"))+
          geom_bar(stat="identity", aes(fill = "P(Target)"))+
          geom_bar(stat="identity")+
          #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          scale_fill_manual("", values = color_vec[1:(length(color_vec)-1)])+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="", y = paste0(""), cex = 1.5) +
          geom_hline(aes(yintercept=ewoc.threshold, linetype ="EWOC"), color = "firebrick", size = 0.8)+
          scale_linetype_manual("", values = "dotted")+
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                size = 2, linetype = "solid"),
                panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                colour = "white"),
                panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                #            colour = "white"),
                plot.margin = unit(c(3, 3, -2, 3), "mm"),
                plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                axis.title = element_text(size=12, face = "bold", family= "sans"),
                axis.text = element_text(size=10, family= "sans"),
                legend.margin = margin(0, 0, -3, 0, "mm"),
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size=12, family= "sans"),
                legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                legend.position = "top")
        #pover
        ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
          geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="Dose Combination", y = paste0("Probability"), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                size = 2, linetype = "solid"),
                panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                colour = "white"),
                panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                colour = "white"),
                plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                axis.title = element_text(size=12, face = "bold", family= "sans"),
                axis.text = element_text(size=10, family= "sans"),
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size=12, family= "sans"),
                legend.text = element_text(size=10, family= "sans"),
                legend.position = "top")

        plist <- list(pover, ptar)


        if(!is.null(file.name) & !is.null(path)){
          if(dir.exists(file.path(path))){
            if(is.null(unit)){
              ggsave(
                filename=paste0(file.name, "_combi.", file.format),
                path=path,
                device=ffint,
                plot=grid.arrange(grobs=plist, nrow=2),
                dpi="retina"
              )
            }else{
              ggsave(
                filename=paste0(file.name, "_combi.", file.format),
                path=path,
                device=ffint,
                plot=grid.arrange(grobs=plist, nrow=2),
                dpi="retina",
                width=width[3],
                height=height[3],
                units = unit
              )
            }
          }
        }
        return(grid.arrange(grobs=plist, nrow=2))

      }

    }else{
      #----------------------------------------------------
      #Combination plot as heatmap
      #----------------------------------------------------
      dose_raw <- rownames(summary)
      ndoses <- length(dose_raw)
      dose1 <- rep(0, times=ndoses)
      dose2 <- rep(0, times=ndoses)
      strdose <- strsplit(dose_raw, split="+", fixed = T)
      for(d in 1:ndoses){
        dose1[d] <- as.numeric(strdose[[d]][1])
        dose2[d] <- as.numeric(strdose[[d]][2])
      }

      lvd1 <- as.numeric(levels(factor(dose1)))
      lvd2 <- as.numeric(levels(factor(dose2)))
      ndose1 <- length(lvd1)
      ndose2 <- length(lvd2)

      datdose1 <- rep(lvd1, each=ndose2)
      datdose2 <- rep(lvd2, times=ndose1)
      name_hmap <- paste0(datdose1,"+", datdose2)

      #Note: not all dose combinations might be in doses of interest!
      #      their tox and ewoc values are left to be NA (later plotted in grey)
      toxs <- vector("numeric", length=ndose1*ndose2)
      names(toxs)<-name_hmap
      ewocs <-vector("numeric", length=ndose1*ndose2)
      names(ewocs)<-name_hmap
      ewodose <- FALSE
      nonewodose <- FALSE
      missdose <- FALSE
      for(i1 in 1:ndose1){
        for(i2 in 1:ndose2){
          dname <- paste0(lvd1[i1], "+", lvd2[i2])
          if(dname %in% dose_raw){
            tox<-summary[dname, 3]
            toxs[(i1-1)*ndose2+i2] <- tox
            if(tox<ewoc.threshold){
              ewocs[(i1-1)*ndose2+i2] <- paste0("P(Over)<", ewoc.threshold)
              ewodose <- TRUE
            }else{
              ewocs[(i1-1)*ndose2+i2] <- paste0("P(Over)\u2265", ewoc.threshold)
              nonewodose <- TRUE
            }
          }else{
            ewocs[(i1-1)*ndose2+i2] <- ""
            toxs[(i1-1)*ndose2+i2] <- ""
            missdose <- TRUE
          }
        }
      }

      #create data frame for plotting
      ldat <- list(
        Dose1=factor(datdose1, levels=unique(datdose1)),
        Dose2=factor(datdose2, levels=unique(datdose2)),
        EWOC=factor(ewocs),
        POD=toxs
      )

#      ldat[["Dose1"]] <-
#      ldat[["Dose2"]] <-
#      ldat[["EWOC"]] <- factor(ewocs)
#      ldat[["POD"]] <- toxs

      datf <- as.data.frame(ldat)

      #color coding depends on whether there are doses not of interest,
      #and whether all doses satisfy or do not satisfy EWOC
      if(missdose&ewodose&nonewodose){
        color_lex <- c("grey95", "forestgreen","firebrick")
        alpha_lex <- c(0, 1, 1)
        names(color_lex)<-c( "", paste0("P(Over)<", ewoc.threshold), paste0("P(Over)\u2265", ewoc.threshold))
        breaks_lex <- c(paste0("P(Over)<", ewoc.threshold), paste0("P(Over)\u2265", ewoc.threshold))
        names(alpha_lex)<-c( "", paste0("P(Over)<", ewoc.threshold), paste0("P(Over)\u2265", ewoc.threshold))
      }else if(!missdose){
        if(ewodose&nonewodose){
          color_lex <- c("forestgreen","firebrick")
          alpha_lex <- c(1, 1)
          names(color_lex)<-c( paste0("P(Over)<", ewoc.threshold), paste0("P(Over)\u2265", ewoc.threshold))
          breaks_lex <- c(paste0("P(Over)<", ewoc.threshold), paste0("P(Over)\u2265", ewoc.threshold))
          names(alpha_lex)<-c( paste0("P(Over)<", ewoc.threshold), paste0("P(Over)\u2265", ewoc.threshold))
        }else if(!ewodose){
          color_lex <- c("firebrick")
          alpha_lex <- c(1)
          names(color_lex)<-c(paste0("P(Over)\u2265", ewoc.threshold))
          breaks_lex <- c(paste0("P(Over)\u2265", ewoc.threshold))
          names(alpha_lex)<-c(paste0("P(Over)\u2265", ewoc.threshold))
        }else if(!nonewodose){
          color_lex <- c("forestgreen")
          alpha_lex <- c(1)
          names(color_lex)<-c( paste0("P(Over)<", ewoc.threshold))
          breaks_lex <- c(paste0("P(Over)<", ewoc.threshold))
          names(alpha_lex)<-c( paste0("P(Over)<", ewoc.threshold))
        }
      }else{
        #there are missing doses
        if(!ewodose){
          color_lex <- c("grey95","firebrick")
          alpha_lex <- c(0,1)
          names(color_lex)<-c("",paste0("P(Over)\u2265", ewoc.threshold))
          breaks_lex <- c(paste0("P(Over)\u2265", ewoc.threshold))
          names(alpha_lex)<-c("",paste0("P(Over)\u2265", ewoc.threshold))
        }else if(!nonewodose){
          color_lex <- c("grey95","forestgreen")
          alpha_lex <- c(0,1)
          names(color_lex)<-c("", paste0("P(Over)<", ewoc.threshold))
          breaks_lex <- c(paste0("P(Over)<", ewoc.threshold))
          names(alpha_lex)<-c("", paste0("P(Over)<", ewoc.threshold))
        }
      }


      heatmap <- ggplot(datf, aes_string(x="Dose1", y="Dose2", fill="EWOC", alpha="EWOC"))+
                  geom_tile(size=1, linetype=1, color="grey95")+
                  scale_fill_manual(breaks=breaks_lex, values = color_lex)+
                  scale_alpha_manual(values = alpha_lex, guide="none")+
                  geom_text( aes_string(label = "POD"), size = 3.5)+
                  labs(#title=
                                 x="Dose 1", y = "Dose 2", cex = 1.5) +
                  theme(
                       panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_blank(),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")

      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=heatmap,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=heatmap,
              dpi="retina",
              width=width[3],
              height=height[3],
              units = unit
            )
          }
        }
      }
      return(heatmap)
    }

  }


}





#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#LOSS ESCALATION
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#----------------------------------------------------
#general plotting function for loss escalation
#----------------------------------------------------

#'@keywords internal
plot_decisions_jointBLRM_loss <- function(
  summary,
  type,
  heatmap,
  loss.weights,
  p.int.probs,
  file.name,
  path,
  file.format,
  width=NULL,
  height=NULL,
  unit=NULL
){
  if(type %in% c("mono1", "mono2", "combi")){
    return(plot_decisions_jointBLRM_losst(summary=summary,
                                          type=type,
                                          heatmap=heatmap,
                                          loss.weights=loss.weights,
                                          p.int.probs = p.int.probs,
                                          file.name=file.name,
                                          path=path,
                                          file.format=file.format,
                                          width=width,
                                          height=height,
                                          unit=unit))
  }else{
    #plot type is "all"

    #Determine if mono1, mono2 and combi doses are in summary,
    #and call the corresponding type-specific plot functions for each
    #detected type.

    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose1 <- rep(0, times=ndoses)
    dose2 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)

    #count number of dose levels per type
    nmon1 <- 0
    nmon2 <- 0
    ncom <- 0
    mono1 <- FALSE
    mono2 <- FALSE
    combi <- FALSE
    for(d in 1:ndoses){
      dose1[d] <- as.numeric(strdose[[d]][1])
      dose2[d] <- as.numeric(strdose[[d]][2])
      if(dose1[d]>0 & dose2[d]==0){
        nmon1 <- nmon1+1
        mono1<-TRUE
      }else if(dose1[d]==0 & dose2[d]>0){
        nmon2 <- nmon2+1
        mono2<-TRUE
      }else if(dose1[d]>0 & dose2[d]>0){
        ncom <- nmon2+1
        combi<-TRUE
      }
    }

    nacc <- mono1+mono2+combi
    plots <- list()

    if(mono1){
      summ <- summary[which(dose1>0 & dose2==0), ]
      plmon1 <- plot_decisions_jointBLRM_losst(summary=summ,
                                               type="mono1",
                                               heatmap = heatmap,
                                               loss.weights=loss.weights,
                                               p.int.probs = p.int.probs,
                                               file.name=file.name,
                                               path=path,
                                               file.format=file.format,
                                               width=width,
                                               height=height,
                                               unit=unit)
      if(nacc==1){
        return(plmon1)
      }else{
        plots[["Mono1"]] <- plmon1
      }
    }

    if(mono2){
      summ <- summary[which(dose2>0 & dose1==0), ]
      plmon2 <- plot_decisions_jointBLRM_losst(summary=summ,
                                               type="mono2",
                                               heatmap = heatmap,
                                               loss.weights=loss.weights,
                                               p.int.probs = p.int.probs,
                                               file.name=file.name,
                                               path=path,
                                               file.format=file.format,
                                               width=width,
                                               height=height,
                                               unit=unit)
      if(nacc==1){
        return(plmon2)
      }else{
        plots[["Mono2"]] <- plmon2
      }
    }

    if(combi){
      summ <- summary[which(dose2>0 & dose1>0), ]
      plcom <- plot_decisions_jointBLRM_losst(summary=summ,
                                              type="combi",
                                              heatmap = heatmap,
                                              loss.weights=loss.weights,
                                              p.int.probs = p.int.probs,
                                              file.name=file.name,
                                              path=path,
                                              file.format=file.format,
                                              width=width,
                                              height=height,
                                              unit=unit)
      if(nacc==1){
        return(plcom)
      }else{
        plots[["Combi"]] <- plcom
      }
    }

    return(plots)

  }

}

#----------------------------------------------------
#type-specific plots for loss escalation
#----------------------------------------------------

#'@keywords internal
plot_decisions_jointBLRM_losst <- function(
  summary,
  type,
  heatmap,
  loss.weights,
  p.int.probs,
  file.name,
  path,
  file.format,
  width=NULL,
  height=NULL,
  unit=NULL
){



  if(file.format=="pdf"){
    ffint <- cairo_pdf
  }else{
    ffint <- file.format
  }

  if(type=="mono1"){
    #----------------------------------------------------
    #Plots for LOSS escalation: Mono 1
    #----------------------------------------------------
    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose1 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)
    for(d in 1:ndoses){
      dose1[d] <- as.numeric(strdose[[d]][1])
    }

    #create data frame for plotting
    ldat <- list()
    loss <- rep(0, times = ndoses)
    vecfill <- rep("P(Over)", times = ndoses)
    for(i in 1:ndoses){
      loss[i] <- sum(loss.weights*summary[i, ])
    }
    ldat[["dose"]] <- factor(dose1, levels=unique(dose1))
    ldat[["under"]] <- summary[,1]
    ldat[["target"]] <- summary[, 2]
    ldat[["over"]] <- summary[, 3]+summary[, 4]
    ldat[["Loss"]]<- loss
    ldat[["labels"]]<- round(loss, digits=5)
    ldat[["vecfill"]]<- vecfill

    datf <- as.data.frame(ldat)

    color_vec <- c( "firebrick", "forestgreen","gold")
    if(p.int.probs){
      punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
        geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill="vecfill"))+
        geom_bar(stat="identity", aes(fill = "P(Under)"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec)+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                       #            colour = "white"),
                       plot.margin = unit(c(3, 3, -2, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.margin = margin(0, 0, -3, 0, "mm"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                       legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")

      plist <- list(pover, ptar, punder)
      p_interval_probs <- grid.arrange(grobs=plist, nrow=3)
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_mono1_int_probs.", file.format),
              path=path,
              device=ffint,
              plot=p_interval_probs,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_mono1_int_probs.", file.format),
              path=path,
              device=ffint,
              plot=p_interval_probs,
              dpi="retina",
              width=width[1],
              height=height[1],
              units = unit
            )

          }
        }
      }
    }

    min.loss <- min(loss)
    ploss <- ggplot(data=datf,aes_string(x="dose", y="Loss", fill="Loss"))+
      geom_bar(stat="identity")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
      scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                          name="Exp. Loss"#, high = "#132B43",low = "#56B1F7"
      )+
      geom_text(aes_string(label = "labels"), y=min.loss*0.6, size=3)+
      labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
        x="Dose", y = paste0("Exp. Loss"), cex = 1.5) +
      geom_hline(aes(yintercept=min.loss, linetype ="Min. Exp. Loss"), color = "steelblue3", size = 0.8)+
      scale_linetype_manual("", values = "dotted")+
      theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                              size = 2, linetype = "solid"),
                     panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                              colour = "white"),
                     panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                              colour = "white"),
                     #plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                     plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                     plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                     axis.title = element_text(size=12, face = "bold", family= "sans"),
                     axis.text = element_text(size=10, family= "sans"),
                     legend.key.size = unit(5, "mm"),
                     legend.title = element_text(size=12, family= "sans"),
                     legend.text = element_text(size=10, family= "sans"),
                     legend.position = "top")
    if(!is.null(file.name) & !is.null(path)){
      if(dir.exists(file.path(path))){
        if(is.null(unit)){
          ggsave(
            filename=paste0(file.name, "_mono1.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina"
          )
        }else{
          ggsave(
            filename=paste0(file.name, "_mono1.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina",
            width=width[1],
            height=height[1],
            units = unit
          )

        }
      }
    }


    if(p.int.probs){
        return(list("ExpLoss"=ploss,
                    "IntervalProb"=p_interval_probs))
    }else{
      return(ploss)
    }
  }else if(type=="mono2"){
    #----------------------------------------------------
    #Plots for LOSS escalation: Mono 2
    #----------------------------------------------------
    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose2 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)
    for(d in 1:ndoses){
      dose2[d] <- as.numeric(strdose[[d]][2])
    }

    #create data frame for plotting
    ldat <- list()
    loss <- rep(0, times = ndoses)
    vecfill <- rep("P(Over)", times = ndoses)
    for(i in 1:ndoses){
      loss[i] <- sum(loss.weights*summary[i, ])
    }
    ldat[["dose"]] <- factor(dose2, levels=unique(dose2))
    ldat[["under"]] <- summary[,1]
    ldat[["target"]] <- summary[, 2]
    ldat[["over"]] <- summary[, 3]+summary[, 4]
    ldat[["Loss"]]<- loss
    ldat[["labels"]]<- round(loss, digits = 5)
    ldat[["vecfill"]] <- vecfill
    datf <- as.data.frame(ldat)

    color_vec <- c( "firebrick", "forestgreen","gold")
    if(p.int.probs){
      punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
        geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill="vecfill"))+
        geom_bar(stat="identity", aes(fill = "P(Under)"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec)+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                       #            colour = "white"),
                       plot.margin = unit(c(3, 3, -2, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.margin = margin(0, 0, -3, 0, "mm"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                       legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")

      plist <- list(pover, ptar, punder)
      p_interval_probs <- grid.arrange(grobs=plist, nrow=3)
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))) {
          if (is.null(unit)) {
            ggsave(
              filename = paste0(file.name, "_mono2_int_probs.", file.format),
              path = path,
              device = ffint,
              plot = p_interval_probs,
              dpi = "retina"
            )
          } else{
            ggsave(
              filename = paste0(file.name, "_mono2_int_probs.", file.format),
              path = path,
              device = ffint,
              plot = p_interval_probs,
              dpi = "retina",
              width = width[2],
              height = height[2],
              units = unit
            )

          }
        }
      }
    }

    min.loss <- min(loss)
    ploss <- ggplot(data=datf,aes_string(x="dose", y="Loss", fill="Loss"))+
      geom_bar(stat="identity")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
      scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                          name="Exp. Loss"#, high = "#132B43",low = "#56B1F7"
      )+
      geom_text(aes_string(label = "labels"), y=min.loss*0.6, size=3)+
      labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
        x="Dose", y = paste0("Exp. Loss"), cex = 1.5) +
      geom_hline(aes(yintercept=min.loss, linetype ="Min. Exp. Loss"), color = "steelblue3", size = 0.8)+
      scale_linetype_manual("", values = "dotted")+
      theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                              size = 2, linetype = "solid"),
                     panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                              colour = "white"),
                     panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                              colour = "white"),
                     #plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                     plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                     plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                     axis.title = element_text(size=12, face = "bold", family= "sans"),
                     axis.text = element_text(size=10, family= "sans"),
                     legend.key.size = unit(5, "mm"),
                     legend.title = element_text(size=12, family= "sans"),
                     legend.text = element_text(size=10, family= "sans"),
                     legend.position = "top")

    if(!is.null(file.name) & !is.null(path)){
      if(dir.exists(file.path(path))){
        if(is.null(unit)){
          ggsave(
            filename=paste0(file.name, "_mono2.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina"
          )
        }else{
          ggsave(
            filename=paste0(file.name, "_mono2.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina",
            width=width[2],
            height=height[2],
            units = unit
          )

        }
      }
    }

    if(p.int.probs){
      return(list("ExpLoss"=ploss,
                  "IntervalProb"=p_interval_probs))
    }else{
      return(ploss)
    }
  }else if(type=="combi"){

    if(!heatmap){
      #----------------------------------------------------
      #Plots for LOSS escalation: Combi
      #----------------------------------------------------
      doses <- rownames(summary)
      ndoses <- length(doses)

      #create data frame for plotting
      ldat <- list()
      loss <- rep(0, times = ndoses)
      vecfill <- rep("P(Over)", times = ndoses)
      for(i in 1:ndoses){
        loss[i] <- sum(loss.weights*summary[i, ])
      }
      ldat[["dose"]] <- factor(doses, levels=unique(doses))
      ldat[["under"]] <- summary[,1]
      ldat[["target"]] <- summary[, 2]
      ldat[["over"]] <- summary[, 3]+summary[, 4]
      ldat[["Loss"]]<- loss
      ldat[["vecfill"]]<- vecfill

      datf <- as.data.frame(ldat)

      color_vec <- c( "firebrick", "forestgreen","gold")
      if(p.int.probs){
        punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
          geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="Dose Combination", y = paste0(""), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                  size = 2, linetype = "solid"),
                         panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                  colour = "white"),
                         panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                  colour = "white"),
                         plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                         plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                         plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                         axis.title = element_text(size=12, face = "bold", family= "sans"),
                         axis.text = element_text(size=10, family= "sans"),
                         legend.key.size = unit(5, "mm"),
                         legend.title = element_text(size=12, family= "sans"),
                         legend.text = element_text(size=10, family= "sans"),
                         legend.position = "top")
        pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill="vecfill"))+
          geom_bar(stat="identity", aes(fill = "P(Under)"))+
          geom_bar(stat="identity", aes(fill = "P(Target)"))+
          geom_bar(stat="identity")+
          #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          scale_fill_manual("", values = color_vec)+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="", y = paste0(""), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                  size = 2, linetype = "solid"),
                         panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                  colour = "white"),
                         panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                         #            colour = "white"),
                         plot.margin = unit(c(3, 3, -2, 3), "mm"),
                         plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                         plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                         axis.title = element_text(size=12, face = "bold", family= "sans"),
                         axis.text = element_text(size=10, family= "sans"),
                         legend.margin = margin(0, 0, -3, 0, "mm"),
                         legend.key.size = unit(5, "mm"),
                         legend.title = element_text(size=12, family= "sans"),
                         legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                         legend.position = "top")
        #pover
        ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
          geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="", y = paste0("Probability"), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                  size = 2, linetype = "solid"),
                         panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                  colour = "white"),
                         panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                  colour = "white"),
                         plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                         plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                         plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                         axis.title = element_text(size=12, face = "bold", family= "sans"),
                         axis.text = element_text(size=10, family= "sans"),
                         legend.key.size = unit(5, "mm"),
                         legend.title = element_text(size=12, family= "sans"),
                         legend.text = element_text(size=10, family= "sans"),
                         legend.position = "top")

        plist <- list(pover, ptar, punder)
        p_interval_probs <- grid.arrange(grobs=plist, nrow=3)
        if(!is.null(file.name) & !is.null(path)){
          if(dir.exists(file.path(path))){
            if(is.null(unit)){
              ggsave(
                filename=paste0(file.name, "_combi_int_probs.", file.format),
                path=path,
                device=ffint,
                plot=p_interval_probs,
                dpi="retina"
              )
            }else{
              ggsave(
                filename=paste0(file.name, "_combi_int_probs.", file.format),
                path=path,
                device=ffint,
                plot=p_interval_probs,
                dpi="retina",
                width=width[3],
                height=height[3],
                units = unit
              )

            }
          }
        }
      }

      min.loss <- min(loss)
      ploss <- ggplot(data=datf,aes_string(x="dose", y="Loss", fill="Loss"))+
        geom_bar(stat="identity")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                            name="Exp. Loss"#, high = "#132B43",low = "#56B1F7"
        )+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0("Exp. Loss"), cex = 1.5) +
        geom_hline(aes(yintercept=min.loss, linetype ="Min. Exp. Loss"), color = "steelblue3", size = 0.8)+
        scale_linetype_manual("", values = "dotted")+
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       #plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=ploss,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=ploss,
              dpi="retina",
              width=width[3],
              height=height[3],
              units = unit
            )

          }
        }
      }

      if(p.int.probs){
        return(list("ExpLoss"=ploss,
                    "IntervalProb"=p_interval_probs))
      }else{
        return(ploss)
      }

    }else{
      #--------------------------------------------------------------------------------------------------
      #Plots for LOSS escalation: Combi as HEATMAP
      #--------------------------------------------------------------------------------------------------
      dose_raw <- rownames(summary)
      ndoses <- length(dose_raw)
      dose1 <- rep(0, times=ndoses)
      dose2 <- rep(0, times=ndoses)
      strdose <- strsplit(dose_raw, split="+", fixed = T)
      for(d in 1:ndoses){
        dose1[d] <- as.numeric(strdose[[d]][1])
        dose2[d] <- as.numeric(strdose[[d]][2])
      }

      lvd1 <- as.numeric(levels(factor(dose1)))
      lvd2 <- as.numeric(levels(factor(dose2)))
      ndose1 <- length(lvd1)
      ndose2 <- length(lvd2)

      datdose1 <- rep(lvd1, each=ndose2)
      datdose2 <- rep(lvd2, times=ndose1)
      name_hmap <- paste0(datdose1,"+", datdose2)

      #Note: not all dose combinations might be in doses of interest!
      #      their tox and ewoc values are left to be NA (later plotted in grey)
      losses <- vector("numeric", length=ndose1*ndose2)
      losslab <- vector("character", length=ndose1*ndose2)
      lossalp <- vector("numeric", length=ndose1*ndose2)
      missdose <- FALSE

      names(losses)<-name_hmap
      for(i1 in 1:ndose1){
        for(i2 in 1:ndose2){
          curr_dname <- paste0(lvd1[i1],"+", lvd2[i2])
          if(curr_dname%in%dose_raw){
            bayesrisk <- sum(loss.weights*summary[curr_dname,])
            losses[(i1-1)*ndose2 + i2] <- bayesrisk
            losslab[(i1-1)*ndose2 + i2] <- paste0(round(bayesrisk, digits=5))
            lossalp[(i1-1)*ndose2 + i2] <-1
          }else{
            missdose <- TRUE
            losses[(i1-1)*ndose2 + i2]<-NA
            losslab[(i1-1)*ndose2 + i2] <- ""
            lossalp[(i1-1)*ndose2 + i2] <-0
          }
        }
      }
      #fill in dummy values for NA doses
      losses[which(is.na(losses))]<-max(losses, na.rm = TRUE)


      #create data frame for plotting
      ldat <- list()
      ldat[["Dose1"]] <-factor(datdose1, levels=unique(datdose1))
      ldat[["Dose2"]] <-factor(datdose2, levels=unique(datdose2))
      ldat[["ExpLoss"]] <- losses
      ldat[["labels"]] <- losslab
      ldat[["alphas"]] <- factor(lossalp)



      datf <- as.data.frame(ldat)
      if(missdose){
        alphalex <- c(0,1)
      }else{
        alphalex <- c(1)
      }


      heatmap <- ggplot(datf, aes_string(x="Dose1", y="Dose2", fill="ExpLoss", alpha="alphas"))+
        geom_tile(size=1, linetype=1, color="grey95")+
        scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                                     name="Expected Loss"#, high = "#132B43",low = "#56B1F7"
                                     )+
        scale_alpha_manual(values = alphalex, guide="none")+
        geom_text( aes_string(label = "labels"), size = 3.5)+
        labs(#title=
          x="Dose 1", y = "Dose 2", cex = 1.5) +
        theme(
          panel.background = element_rect(fill = "grey95", colour = "white",
                                                   size = 2, linetype = "solid"),
          panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                   colour = "white"),
          panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                   colour = "white"),
          plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
          plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
          axis.title = element_text(size=12, face = "bold", family= "sans"),
          axis.text = element_text(size=10, family= "sans"),
          legend.key.size = unit(5, "mm"),
          legend.title = element_text(size=12, family= "sans"),
          legend.text = element_text(size=10, family= "sans"),
          legend.position = "top")
      #heatmap
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=heatmap,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=heatmap,
              dpi="retina",
              width=width[3],
              height=height[3],
              units = unit
            )

          }
        }
      }

      return(heatmap)

    }


  }

}


#----------------------------------------------------
#----------------------------------------------------
#Dynamic LOSS ESCALATION
#----------------------------------------------------
#----------------------------------------------------

#----------------------------------------------------
#general plotting fun. for dynamic loss escalation
#----------------------------------------------------

#'@keywords internal
plot_decisions_jointBLRM_dloss <- function(
  summary,
  type,
  probs.ref,
  heatmap,
  dloss.weights,
  p.int.probs,
  file.name,
  path,
  file.format,
  width=NULL,
  height=NULL,
  unit=NULL
){
  if(type %in% c("mono1", "mono2", "combi")){
    return(plot_decisions_jointBLRM_dlosst(summary=summary,
                                          type=type,
                                          probs.ref=probs.ref,
                                          heatmap=heatmap,
                                          dloss.weights=dloss.weights,
                                          p.int.probs = p.int.probs,
                                          file.name=file.name,
                                          path=path,
                                          file.format=file.format,
                                          width=width,
                                          height=height,
                                          unit=unit))
  }else{
    #plot type is "all"

    #Determine if mono1, mono2 and combi doses are in summary,
    #and call the corresponding type-specific plot functions for each
    #detected type.

    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose1 <- rep(0, times=ndoses)
    dose2 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)

    #count number of dose levels per type
    nmon1 <- 0
    nmon2 <- 0
    ncom <- 0
    mono1 <- FALSE
    mono2 <- FALSE
    combi <- FALSE
    for(d in 1:ndoses){
      dose1[d] <- as.numeric(strdose[[d]][1])
      dose2[d] <- as.numeric(strdose[[d]][2])
      if(dose1[d]>0 & dose2[d]==0){
        nmon1 <- nmon1+1
        mono1<-TRUE
      }else if(dose1[d]==0 & dose2[d]>0){
        nmon2 <- nmon2+1
        mono2<-TRUE
      }else if(dose1[d]>0 & dose2[d]>0){
        ncom <- nmon2+1
        combi<-TRUE
      }
    }

    nacc <- mono1+mono2+combi
    plots <- list()

    if(mono1){
      summ <- summary[which(dose1>0 & dose2==0), ]
      plmon1 <- plot_decisions_jointBLRM_dlosst(summary=summ,
                                               type="mono1",
                                               heatmap = heatmap,
                                               probs.ref=probs.ref,
                                               dloss.weights=dloss.weights,
                                               p.int.probs = p.int.probs,
                                               file.name=file.name,
                                               path=path,
                                               file.format=file.format,
                                               width=width,
                                               height=height,
                                               unit=unit)
      if(nacc==1){
        return(plmon1)
      }else{
        plots[["Mono1"]] <- plmon1
      }
    }

    if(mono2){
      summ <- summary[which(dose2>0 & dose1==0), ]
      plmon2 <- plot_decisions_jointBLRM_dlosst(summary=summ,
                                               type="mono2",
                                               heatmap = heatmap,
                                               probs.ref=probs.ref,
                                               dloss.weights=dloss.weights,
                                               p.int.probs = p.int.probs,
                                               file.name=file.name,
                                               path=path,
                                               file.format=file.format,
                                               width=width,
                                               height=height,
                                               unit=unit)
      if(nacc==1){
        return(plmon2)
      }else{
        plots[["Mono2"]] <- plmon2
      }
    }

    if(combi){
      summ <- summary[which(dose2>0 & dose1>0), ]
      plcom <- plot_decisions_jointBLRM_dlosst(summary=summ,
                                              type="combi",
                                              heatmap = heatmap,
                                              probs.ref=probs.ref,
                                              dloss.weights=dloss.weights,
                                              p.int.probs = p.int.probs,
                                              file.name=file.name,
                                              path=path,
                                              file.format=file.format,
                                              width=width,
                                              height=height,
                                              unit=unit)
      if(nacc==1){
        return(plcom)
      }else{
        plots[["Combi"]] <- plcom
      }
    }

    return(plots)

  }
}


#----------------------------------------------------
#type-specific plots for dynamic loss escalation
#----------------------------------------------------

#'@keywords internal
plot_decisions_jointBLRM_dlosst <- function(
  summary,
  type,
  heatmap,
  probs.ref,
  dloss.weights,
  p.int.probs,
  file.name,
  path,
  file.format,
  width=NULL,
  height=NULL,
  unit=NULL
){


  if(file.format=="pdf"){
    ffint <- cairo_pdf
  }else{
    ffint <- file.format
  }

  if(type=="mono1"){
    #--------------------------------------------------------------------------------------------------
    #Dynamic LOSS: Mono1
    #--------------------------------------------------------------------------------------------------
    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose1 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)
    for(d in 1:ndoses){
      dose1[d] <- as.numeric(strdose[[d]][1])
    }

    #compute dynamic loss weights
    dlw <- dloss.weights[1,]*probs.ref[1,1] + dloss.weights[2,]*probs.ref[1,2] +
      dloss.weights[3,]*probs.ref[1,3] + dloss.weights[4,]*probs.ref[1,4]
    display.dlw <- round(dlw, digits=3)

    #create data frame for plotting
    ldat <- list()
    loss <- rep(0, times = ndoses)
    vecfill <- rep("P(Over)", times = ndoses)
    for(i in 1:ndoses){
      loss[i] <- sum(dlw*summary[i, ])
    }
    ldat[["dose"]] <- factor(dose1, levels=unique(dose1))
    ldat[["under"]] <- summary[,1]
    ldat[["target"]] <- summary[, 2]
    ldat[["over"]] <- summary[, 3]+summary[, 4]
    ldat[["Loss"]]<- loss
    ldat[["labels"]]<- round(loss, digits = 5)
    ldat[["vecfill"]]<- vecfill

    datf <- as.data.frame(ldat)

    color_vec <- c( "firebrick", "forestgreen","gold")
    if(p.int.probs){
      punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
        geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill="vecfill"))+
        geom_bar(stat="identity", aes(fill = "P(Under)"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec)+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                       #            colour = "white"),
                       plot.margin = unit(c(3, 3, -2, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.margin = margin(0, 0, -3, 0, "mm"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                       legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")

      plist <- list(pover, ptar, punder)
      p_interval_probs <- grid.arrange(grobs=plist, nrow=3)
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_mono1_int_probs.", file.format),
              path=path,
              device=ffint,
              plot=p_interval_probs,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_mono1_int_probs.", file.format),
              path=path,
              device=ffint,
              plot=p_interval_probs,
              dpi="retina",
              width=width[1],
              height=height[1],
              units = unit
            )

          }
        }
      }
    }

    min.loss <- min(loss)
    ploss <- ggplot(data=datf,aes_string(x="dose", y="Loss", fill="Loss"))+
      geom_bar(stat="identity")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
      scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                          name="Exp. Loss"#, high = "#132B43",low = "#56B1F7"
      )+
      geom_text(aes_string(label = "labels"), y=min.loss*0.6, size=3)+
      labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
        subtitle=paste0("Dynamic Loss Weights: (", display.dlw[1], ", ", display.dlw[2], ", ",
                        display.dlw[3], ", ", display.dlw[4], ")"),
        x="Dose", y = paste0("Exp. Loss"), cex = 1.5) +
      geom_hline(aes(yintercept=min.loss, linetype ="Min. Exp. Loss"), color = "steelblue3", size = 0.8)+
      scale_linetype_manual("", values = "dotted")+
      theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                              size = 2, linetype = "solid"),
                     panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                              colour = "white"),
                     panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                              colour = "white"),
                     #plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                     plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                     plot.subtitle = element_text(hjust = 0.5, size = 14, family= "sans"),
                     axis.title = element_text(size=12, face = "bold", family= "sans"),
                     axis.text = element_text(size=10, family= "sans"),
                     legend.key.size = unit(5, "mm"),
                     legend.title = element_text(size=12, family= "sans"),
                     legend.text = element_text(size=10, family= "sans"),
                     legend.position = "top")
    if(!is.null(file.name) & !is.null(path)){
      if(dir.exists(file.path(path))){
        if(is.null(unit)){
          ggsave(
            filename=paste0(file.name, "_mono1.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina"
          )
        }else{
          ggsave(
            filename=paste0(file.name, "_mono1.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina",
            width=width[1],
            height=height[1],
            units = unit
          )

        }
      }
    }

    if(p.int.probs){
      return(list("ExpLoss"=ploss,
                  "IntervalProb"=p_interval_probs))
    }else{
      return(ploss)
    }


  }else if(type=="mono2"){
    #--------------------------------------------------------------------------------------------------
    #Dynamic LOSS: Mono2
    #--------------------------------------------------------------------------------------------------
    dose_raw <- rownames(summary)
    ndoses <- length(dose_raw)
    dose2 <- rep(0, times=ndoses)
    strdose <- strsplit(dose_raw, split="+", fixed = T)
    for(d in 1:ndoses){
      dose2[d] <- as.numeric(strdose[[d]][2])
    }

    #compute dynamic loss weights
    dlw <- dloss.weights[1,]*probs.ref[2,1] + dloss.weights[2,]*probs.ref[2,2] +
      dloss.weights[3,]*probs.ref[2,3] + dloss.weights[4,]*probs.ref[2,4]
    display.dlw <- round(dlw, digits=3)

    #create data frame for plotting
    ldat <- list()
    loss <- rep(0, times = ndoses)
    vecfill <- rep("P(Over)", times = ndoses)
    for(i in 1:ndoses){
      loss[i] <- sum(dlw*summary[i, ])
    }
    ldat[["dose"]] <- factor(dose2, levels=unique(dose2))
    ldat[["under"]] <- summary[,1]
    ldat[["target"]] <- summary[, 2]
    ldat[["over"]] <- summary[, 3]+summary[, 4]
    ldat[["Loss"]]<- loss
    ldat[["labels"]]<- round(loss, digits = 5)

    ldat[["vecfill"]]<- vecfill
    datf <- as.data.frame(ldat)

    color_vec <- c( "firebrick", "forestgreen","gold")
    if(p.int.probs){
      punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
        geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")
      pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill="vecfill"))+
        geom_bar(stat="identity", aes(fill = "P(Under)"))+
        geom_bar(stat="identity", aes(fill = "P(Target)"))+
        geom_bar(stat="identity")+
        #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_manual("", values = color_vec)+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0(""), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                       #            colour = "white"),
                       plot.margin = unit(c(3, 3, -2, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.margin = margin(0, 0, -3, 0, "mm"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                       legend.position = "top")
      #pover
      ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
        geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
          x="", y = paste0("Probability"), cex = 1.5) +
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")

      plist <- list(pover, ptar, punder)
      p_interval_probs <- grid.arrange(grobs=plist, nrow=3)
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_mono2_int_probs.", file.format),
              path=path,
              device=ffint,
              plot=p_interval_probs,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_mono2_int_probs.", file.format),
              path=path,
              device=ffint,
              plot=p_interval_probs,
              dpi="retina",
              width=width[2],
              height=height[2],
              units = unit
            )

          }
        }
      }
    }

    min.loss <- min(loss)
    ploss <- ggplot(data=datf,aes_string(x="dose", y="Loss", fill="Loss"))+
      geom_bar(stat="identity")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
      scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                          name="Exp. Loss"#, high = "#132B43",low = "#56B1F7"
      )+
      geom_text(aes_string(label = "labels"), y=min.loss*0.6, size=3)+
      labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
        subtitle=paste0("Dynamic Loss Weights: (", display.dlw[1], ", ", display.dlw[2], ", ",
                        display.dlw[3], ", ", display.dlw[4], ")"),
        x="Dose", y = paste0("Exp. Loss"), cex = 1.5) +
      geom_hline(aes(yintercept=min.loss, linetype ="Min. Exp. Loss"), color = "steelblue3", size = 0.8)+
      scale_linetype_manual("", values = "dotted")+
      theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                              size = 2, linetype = "solid"),
                     panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                              colour = "white"),
                     panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                              colour = "white"),
                     #plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                     plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                     plot.subtitle = element_text(hjust = 0.5, size = 14, family= "sans"),
                     axis.title = element_text(size=12, face = "bold", family= "sans"),
                     axis.text = element_text(size=10, family= "sans"),
                     legend.key.size = unit(5, "mm"),
                     legend.title = element_text(size=12, family= "sans"),
                     legend.text = element_text(size=10, family= "sans"),
                     legend.position = "top")
    if(!is.null(file.name) & !is.null(path)){
      if(dir.exists(file.path(path))){
        if(is.null(unit)){
          ggsave(
            filename=paste0(file.name, "_mono2.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina"
          )
        }else{
          ggsave(
            filename=paste0(file.name, "_mono2.", file.format),
            path=path,
            device=ffint,
            plot=ploss,
            dpi="retina",
            width=width[2],
            height=height[2],
            units = unit
          )

        }
      }
    }

    if(p.int.probs){
      return(list("ExpLoss"=ploss,
                  "IntervalProb"=p_interval_probs))
    }else{
      return(ploss)
    }


  }else if(type=="combi"){
    if(!heatmap){
      #--------------------------------------------------------------------------------------------------
      #Dynamic LOSS: Combi, no heatmap
      #--------------------------------------------------------------------------------------------------
      doses <- rownames(summary)
      ndoses <- length(doses)

      #create data frame for plotting
      ldat <- list()
      loss <- rep(0, times = ndoses)
      vecfill <- rep("P(Over)", times = ndoses)
      #compute dynamic loss weights
      dlw <- dloss.weights[1,]*probs.ref[3,1] + dloss.weights[2,]*probs.ref[3,2] +
        dloss.weights[3,]*probs.ref[3,3] + dloss.weights[4,]*probs.ref[3,4]
      display.dlw <- round(dlw, digits=3)
      for(i in 1:ndoses){
        loss[i] <- sum(dlw*summary[i, ])
      }
      ldat[["dose"]] <- factor(doses, levels=unique(doses))
      ldat[["under"]] <- summary[,1]
      ldat[["target"]] <- summary[, 2]
      ldat[["over"]] <- summary[, 3]+summary[, 4]
      ldat[["Loss"]]<- loss
      ldat[["vecfill"]]<- vecfill
      datf <- as.data.frame(ldat)

      color_vec <- c( "firebrick", "forestgreen","gold")
      if(p.int.probs){
        punder <- ggplot(data=datf,aes_string(x="dose", y="under"))+
          geom_bar(stat="identity", fill = "gold")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="Dose Combination", y = paste0(""), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                  size = 2, linetype = "solid"),
                         panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                  colour = "white"),
                         panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                  colour = "white"),
                         plot.margin = unit(c(-1, 3, 3, 3), "mm"),
                         plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                         plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                         axis.title = element_text(size=12, face = "bold", family= "sans"),
                         axis.text = element_text(size=10, family= "sans"),
                         legend.key.size = unit(5, "mm"),
                         legend.title = element_text(size=12, family= "sans"),
                         legend.text = element_text(size=10, family= "sans"),
                         legend.position = "top")
        pover <- ggplot(data=datf,aes_string(x="dose", y="over", fill="vecfill"))+
          geom_bar(stat="identity", aes(fill = "P(Under)"))+
          geom_bar(stat="identity", aes(fill = "P(Target)"))+
          geom_bar(stat="identity")+
          #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          scale_fill_manual("", values = color_vec)+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="", y = paste0(""), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                  size = 2, linetype = "solid"),
                         panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                  colour = "white"),
                         panel.grid.minor = element_blank(),#element_line(size = 0.2, linetype = 'solid',
                         #            colour = "white"),
                         plot.margin = unit(c(3, 3, -2, 3), "mm"),
                         plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                         plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                         axis.title = element_text(size=12, face = "bold", family= "sans"),
                         axis.text = element_text(size=10, family= "sans"),
                         legend.margin = margin(0, 0, -3, 0, "mm"),
                         legend.key.size = unit(5, "mm"),
                         legend.title = element_text(size=12, family= "sans"),
                         legend.text = element_text(size=10, family= "sans", margin = margin(0, 4, 0, -4.5, "pt")),
                         legend.position = "top")
        #pover
        ptar <- ggplot(data=datf,aes_string(x="dose", y="target"))+
          geom_bar(stat="identity", fill = "forestgreen")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
          labs(#title=paste0("Medians and credible intervals of the prior toxicities"),
            x="", y = paste0("Probability"), cex = 1.5) +
          theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                  size = 2, linetype = "solid"),
                         panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                  colour = "white"),
                         panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                  colour = "white"),
                         plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                         plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                         plot.subtitle = element_text(hjust = 0.5, size = 18, family= "sans"),
                         axis.title = element_text(size=12, face = "bold", family= "sans"),
                         axis.text = element_text(size=10, family= "sans"),
                         legend.key.size = unit(5, "mm"),
                         legend.title = element_text(size=12, family= "sans"),
                         legend.text = element_text(size=10, family= "sans"),
                         legend.position = "top")

        plist <- list(pover, ptar, punder)
        p_interval_probs <- grid.arrange(grobs=plist, nrow=3)
        if(!is.null(file.name) & !is.null(path)){
          if(dir.exists(file.path(path))){
            if(is.null(unit)){
              ggsave(
                filename=paste0(file.name, "_combi_int_probs.", file.format),
                path=path,
                device=ffint,
                plot=p_interval_probs,
                dpi="retina"
              )
            }else{
              ggsave(
                filename=paste0(file.name, "_combi_int_probs.", file.format),
                path=path,
                device=ffint,
                plot=p_interval_probs,
                dpi="retina",
                width=width[3],
                height=height[3],
                units = unit
              )

            }
          }
        }
      }

      min.loss <- min(loss)
      ploss <- ggplot(data=datf,aes_string(x="dose", y="Loss", fill="Loss"))+
        geom_bar(stat="identity")+ #scale_x_continuous(breaks = data$doses1, trans = "log2")+
        scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                            name="Exp. Loss"#, high = "#132B43",low = "#56B1F7"
        )+
        labs(subtitle=paste0("Dynamic Loss Weights: (", display.dlw[1], ", ", display.dlw[2], ", ",
                                      display.dlw[3], ", ", display.dlw[4], ")"),
                      #title=paste0("Medians and credible intervals of the prior toxicities"),
          x="Dose", y = paste0("Exp. Loss"), cex = 1.5) +
        geom_hline(aes(yintercept=min.loss, linetype ="Min. Exp. Loss"), color = "steelblue3", size = 0.8)+
        scale_linetype_manual("", values = "dotted")+
        theme(panel.background = element_rect(fill = "grey95", colour = "white",
                                                                size = 2, linetype = "solid"),
                       panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                                colour = "white"),
                       panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                                colour = "white"),
                       #plot.margin = unit(c(-1, 3, -1, 3), "mm"),
                       plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
                       plot.subtitle = element_text(hjust = 0.5, size = 14, family= "sans"),
                       axis.title = element_text(size=12, face = "bold", family= "sans"),
                       axis.text = element_text(size=10, family= "sans"),
                       legend.key.size = unit(5, "mm"),
                       legend.title = element_text(size=12, family= "sans"),
                       legend.text = element_text(size=10, family= "sans"),
                       legend.position = "top")
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=ploss,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=ploss,
              dpi="retina",
              width=width[3],
              height=height[3],
              units = unit
            )

          }
        }
      }

      if(p.int.probs){
        return(list("ExpLoss"=ploss,
                    "IntervalProb"=p_interval_probs))
      }else{
        return(ploss)
      }


    }else{
      #--------------------------------------------------------------------------------------------------
      #Dynamic LOSS: combi with HEATMAP
      #--------------------------------------------------------------------------------------------------
      dose_raw <- rownames(summary)
      ndoses <- length(dose_raw)
      dose1 <- rep(0, times=ndoses)
      dose2 <- rep(0, times=ndoses)
      strdose <- strsplit(dose_raw, split="+", fixed = TRUE)
      for(d in 1:ndoses){
        dose1[d] <- as.numeric(strdose[[d]][1])
        dose2[d] <- as.numeric(strdose[[d]][2])
      }

      lvd1 <- as.numeric(levels(factor(dose1)))
      lvd2 <- as.numeric(levels(factor(dose2)))
      ndose1 <- length(lvd1)
      ndose2 <- length(lvd2)

      datdose1 <- rep(lvd1, each=ndose2)
      datdose2 <- rep(lvd2, times=ndose1)
      name_hmap <- paste0(datdose1,"+", datdose2)

      #Note: not all dose combinations might be in doses of interest!
      #      their tox and ewoc values are left to be NA (later plotted in grey)
      losses <- vector("numeric", length=ndose1*ndose2)
      losslab <- vector("character", length=ndose1*ndose2)
      lossalp <- vector("numeric", length=ndose1*ndose2)
      missdose <- FALSE

      #compute dynamic loss weights
      dlw <- dloss.weights[1,]*probs.ref[3,1] + dloss.weights[2,]*probs.ref[3,2] +
        dloss.weights[3,]*probs.ref[3,3] + dloss.weights[4,]*probs.ref[3,4]
      display.dlw <- round(dlw, digits=3)

      names(losses)<-name_hmap
      for(i1 in 1:ndose1){
        for(i2 in 1:ndose2){
          curr_dname <- paste0(lvd1[i1],"+", lvd2[i2])
          if(curr_dname%in%dose_raw){
            bayesrisk <- sum(dlw*summary[curr_dname,])
            losses[(i1-1)*ndose2 + i2] <- bayesrisk
            losslab[(i1-1)*ndose2 + i2] <- paste0(round(bayesrisk, digits=5))
            lossalp[(i1-1)*ndose2 + i2] <-1
          }else{
            missdose <- TRUE
            losses[(i1-1)*ndose2 + i2]<-NA
            losslab[(i1-1)*ndose2 + i2] <- ""
            lossalp[(i1-1)*ndose2 + i2] <-0
          }
        }
      }
      #fill in dummy values for NA doses
      losses[which(is.na(losses))]<-max(losses, na.rm = TRUE)


      #create data frame for plotting
      ldat <- list()
      ldat[["Dose1"]] <-factor(datdose1, levels=unique(datdose1))
      ldat[["Dose2"]] <-factor(datdose2, levels=unique(datdose2))
      ldat[["ExpLoss"]] <- losses
      ldat[["labels"]] <- losslab
      ldat[["alphas"]] <- factor(lossalp)



      datf <- as.data.frame(ldat)
      if(missdose){
        alphalex <- c(0,1)
      }else{
        alphalex <- c(1)
      }


      heatmap <- ggplot(datf, aes_string(x="Dose1", y="Dose2", fill="ExpLoss", alpha="alphas"))+
        geom_tile(size=1, linetype=1, color="grey95")+
        scale_fill_gradient(low="powderblue", high="steelblue4", breaks=breaks_extended(3),
                                     name="Expected Loss"#, high = "#132B43",low = "#56B1F7"
        )+
        scale_alpha_manual(values = alphalex, guide="none")+
        geom_text( aes_string(label = "labels"), size = 3.5)+
        labs(subtitle=paste0("Dynamic Loss Weights: (", display.dlw[1], ", ", display.dlw[2], ", ",
                                      display.dlw[3], ", ", display.dlw[4], ")"),
          x="Dose 1", y = "Dose 2", cex = 1.5) +
        theme(
          panel.background = element_rect(fill = "grey95", colour = "white",
                                                   size = 2, linetype = "solid"),
          panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                                   colour = "white"),
          panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                   colour = "white"),
          plot.title = element_text(hjust = 0.5, size = 20, family= "sans"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, family= "sans"),
          axis.title = element_text(size=12, face = "bold", family= "sans"),
          axis.text = element_text(size=10, family= "sans"),
          legend.key.size = unit(5, "mm"),
          legend.title = element_text(size=12, family= "sans"),
          legend.text = element_text(size=10, family= "sans"),
          legend.position = "top")
      #heatmap
      if(!is.null(file.name) & !is.null(path)){
        if(dir.exists(file.path(path))){
          if(is.null(unit)){
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=heatmap,
              dpi="retina"
            )
          }else{
            ggsave(
              filename=paste0(file.name, "_combi.", file.format),
              path=path,
              device=ffint,
              plot=heatmap,
              dpi="retina",
              width=width[3],
              height=height[3],
              units = unit
            )

          }
        }
      }

      return(heatmap)


    }


  }

}
